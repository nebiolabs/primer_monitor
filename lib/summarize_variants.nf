nextflow.enable.dsl=2

ref = params.ref 
ref = file(ref).toAbsolutePath()
params.prev_json=

prev_json = file(params.prev_json, checkIfExists: true).toAbsolutePath()

params.pangolin_path=

pangolin_path = file(params.pangolin_path).toAbsolutePath()

ncov_path = '/mnt/home/mcampbell/src/ncov-ingest'
primer_monitor_path = '/mnt/bioinfo/prg/primer_monitor'
output_path = '/mnt/hpc_scratch/primer_monitor'
pangolin_path = 

process download_data {
    // Downloads the full dataset
    cpus 16
    penv 'smp'
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
    penv 'smp'
    conda "python=3.9 zstd"

    output:
    file '*.json'

    shell:
    '''
    date_today=$(date +%Y-%m-%d)

    python3 !{primer_monitor_path}/lib/filter_duplicates.py <(zstd -d --long=30 < !{prev_json}) <(zstd -d --long=30 < !{full_json}) > ${date_today}.json

    find !{output_path} -maxdepth 1 -mtime +5 -type f -name "*.full_json*"  -delete
    '''
 
}

process transform_data {
    cpus 1
    penv 'smp'
    conda "python=3.9 regex fsspec pandas typing"

    input:
    file gisaid_json

    output:
    tuple file('*.metadata'), file('*.fasta')

    shell:
    '''
    date_today=$(date +%Y-%m-%d)

    !{ncov_path}/bin/transform-gisaid --output-metadata ${date_today}.metadata --output-fasta ${date_today}.fasta --output-additional-info ${date_today}.info ${date_today}.*.json 
    '''

}

process align {
    cpus 16
    penv 'smp'
    conda "minimap2=2.17 sed python=3.9 samtools=1.11"
    publishDir "${output_path}", mode: 'copy', pattern: '*.bam', overwrite: true

    input:
    tuple file(metadata), file(fasta)

    output:
    tuple file('*.metadata'), file('*.tsv')

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
    penv 'smp'
    publishDir "${output_path}", mode: 'copy'
    errorStrategy 'retry' 
    maxRetries 10
    maxForks 1
    input:
        tuple file(metadata), file(tsv)
    output:
        file '*.complete'
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
    cpus 8
    penv 'smp'
    conda "pangolin"
    input:
        tuple file(metadata), file(fasta)
    output:
        file "*.csv"
    shell:
    '''
    !{primer_monitor_path}/lib/pangolin_calls/run_pangolin.sh !{fasta} 8
    '''
}

process load_pangolin_data {
    cpus 1
    penv 'smp'
    input:
        file csv
        file complete
        //the .complete is only here to make sure this happens *after* the main DB load
    output:
        file '*.complete'
    shell:
    '''
    !{primer_monitor_path}/lib/pangolin_calls/update_fasta_records.sh !{csv}
    '''
}

process recalculate_database_views {
    cpus 1
    penv 'smp'
    publishDir "${output_path}", mode: 'copy'
    errorStrategy 'retry' 
    maxRetries 2
    input:
        file everything
        file everything_pangolin
    shell:
    '''
    # recalculate all the views at the end to save time
    RAILS_ENV=production ruby /mnt/bioinfo/prg/primer_monitor/upload.rb --skip_data_import && touch refresh_complete.txt
    '''
}

workflow {
    extract_new_records()
    transform_data(extract_new_records.out.splitText(file: true, by: 10000))
    align(transform_data.out)
    load_to_db(align.out)
    pangolin_calls(transform_data.out)
    load_pangolin_data(pangolin_calls.out, load_to_db.out)
    recalculate_database_views(load_to_db.out.collect(), load_pangolin_data.out.collect())
}
