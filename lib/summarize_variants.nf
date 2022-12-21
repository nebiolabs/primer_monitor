ref = params.ref 
ref = file(ref).toAbsolutePath()
params.prev_json=

prev_json = file(params.prev_json, checkIfExists: true).toAbsolutePath()

params.pangolin_path=

pangolin_path = file(params.pangolin_path).toAbsolutePath()

threads = 8


ncov_path = '/mnt/home/mcampbell/src/ncov-ingest'
primer_monitor_path = '/mnt/bioinfo/prg/primer_monitor'
output_path = '/mnt/hpc_scratch/primer_monitor'
pangolin_path = 

process download_data {
    // Downloads the full dataset
    cpus 16
    conda "curl xz zstd"
    errorStrategy 'retry' 
    maxRetries 2
    publishDir "${output_path}", mode: 'link', pattern: '*.full_json.zst', overwrite: true
    // mode "link" assumes that the output path is on the same disk as the work directory, switch to copy if not

    output:
        file('*.full_json.zst') into downloaded_data

    shell:
    '''
    date_today=$(date +%Y-%m-%d)
    source !{primer_monitor_path}/.env
    curl -u $USER:$PASSWORD $URL > tmp.json.xz
    xz -d < tmp.json.xz | zstd --long=30 --ultra -22 -T!{task.cpus} > ${date_today}.full_json.zst
    rm tmp.json.xz
    '''

}

process extract_new_records {
    // Keeps only new records added since previous run
    cpus 1
    conda "python=3.9 zstd"

    input:
        file(full_json) from downloaded_data
    output:
        file('*.json') into filtered_data


    shell:
    '''
    date_today=$(date +%Y-%m-%d)

    python3 !{primer_monitor_path}/lib/filter_duplicates.py <(zstd -d --long=30 < !{prev_json}) <(zstd -d --long=30 < !{full_json}) > ${date_today}.json

    find !{output_path} -maxdepth 1 -mtime +5 -type f -name "*.full_json*"  -delete
    '''
 
}

process transform_data {
    cpus 1
    conda "regex fsspec pandas typing"

    input:
        file(gisaid_json) from filtered_data.splitText(file: true, by: 10000)
    output:
        tuple file('*.metadata'), file('*.fasta') into transformed_data
        file('*.fasta') into transformed_data_for_pangolin


    shell:
    '''
    date_today=$(date +%Y-%m-%d)

    !{ncov_path}/bin/transform-gisaid --output-metadata ${date_today}.metadata --output-fasta ${date_today}.fasta --output-additional-info ${date_today}.info ${date_today}.*.json 
    '''

}

process align {
    cpus 16
    conda "minimap2=2.17 sed python=3.9 samtools=1.11"
    publishDir "${output_path}", mode: 'copy', pattern: '*.bam', overwrite: true

    input:
        tuple file(metadata), file(fasta) from transformed_data
    output:
        tuple file('*.metadata'), file('*.tsv') into metadata_plus_variants

    shell:
    '''
        date_today=$(date +%Y-%m-%d)

        sed -E 's/ /_/g' !{fasta} \
        | minimap2 -r 10000 --score-N=0 -t !{task.cpus} --eqx -x map-ont -a !{ref} /dev/stdin \
        | samtools view -h -F 2304 /dev/stdin \
        | tee >(samtools view -b -o ${date_today}.$(basename $PWD).bam /dev/stdin) \
        | python3 !{primer_monitor_path}/lib/parse_alignments.py > $(basename $PWD).tsv

        mv !{metadata} $(basename $PWD).metadata

    '''
}

process load_to_db {
    cpus 1
    publishDir "${output_path}", mode: 'copy'
    errorStrategy 'retry' 
    maxRetries 10
    maxForks 1
    input:
        tuple file(metadata), file(tsv) from metadata_plus_variants
    output:
        file('*.complete') into complete_metadata_files
        file('*.complete') into complete_metadata_files_pangolin
    shell:
    '''
    RAILS_ENV=production ruby /mnt/bioinfo/prg/primer_monitor/upload.rb \
        --skip_view_rebuild \
        --metadata_tsv !{metadata} \
        --variants_tsv !{tsv} \
        && mv !{metadata} !{metadata}.complete
    '''
}

process pangolin_calls {
    cpus 1
    //the actual pangolin runs are qsubbed separately, this 1 cpu is for the main process
    input:
        file(fasta) from transformed_data_for_pangolin
    output:
        file("*.csv") into pangolin_lineage_data
    shell:
    '''
    !{primer_monitor_path}/lib/pangolin_calls/qsub_pangolin.sh !{fasta} !{pangolin_path} !{threads}
    '''
}

process load_pangolin_data {
    cpus 1
    input:
        file(csv) from pangolin_lineage_data
        file(complete) from complete_metadata_files_pangolin
        //the .complete is only here to make sure this happens *after* the main DB load
    output:
        file('*.complete') into complete_files_pangolin
    shell:
    '''
    !{primer_monitor_path}/lib/pangolin_calls/update_fasta_records.sh !{csv}
    '''
}

process recalculate_database_views {
    cpus 1
    publishDir "${output_path}", mode: 'copy'
    errorStrategy 'retry' 
    maxRetries 2
    input:
        file(everything) from complete_metadata_files.collect()
        file(everything_pangolin) from complete_files_pangolin.collect()
    shell:
    '''
    # recalculate all the views at the end to save time
    RAILS_ENV=production ruby /mnt/bioinfo/prg/primer_monitor/upload.rb --skip_data_import && touch refresh_complete.txt
    '''
}
