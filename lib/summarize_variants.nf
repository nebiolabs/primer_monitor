nextflow.enable.dsl = 2

ref = params.ref
ref = file(ref).toAbsolutePath()
params.prev_json=

params.flag_path='/mnt/hpc_scratch/primer_monitor'

prev_json = file(params.prev_json, checkIfExists: true).toAbsolutePath()

ncov_path = '/mnt/home/mcampbell/src/ncov-ingest'
primer_monitor_path = '/mnt/bioinfo/prg/primer_monitor'
output_path = '/mnt/hpc_scratch/primer_monitor'

params.pangolin_version_path =
params.pangolin_data_version_path =

pangolin_version = file(params.pangolin_version_path).text
pangolin_data_version = file(params.pangolin_data_version_path).text

process download_data {
    // Downloads the full dataset
    cpus 16
    penv 'smp'
    conda "ncbi-datasets-cli unzip zstd"
    errorStrategy 'retry'
    maxRetries 2
    publishDir "${output_path}", mode: 'link', pattern: '*.zst', overwrite: true
    // mode "link" assumes that the output path is on the same disk as the work directory, switch to copy if not

    output:
        tuple file('*.metadata.zst'), file('*.sequences.zst')

    shell:
    '''
    date_today=$(date +%Y-%m-%d)
    source !{primer_monitor_path}/.env
    datasets download virus genome taxon SARS-CoV-2 --complete-only --host human --filename tmp.zip
    unzip tmp.zip
    zstd --long=30 --ultra -22 -T!{task.cpus} ncbi_dataset/data/data_report.jsonl -o ${date_today}.metadata.zst
    zstd --long=30 --ultra -22 -T!{task.cpus} ncbi_dataset/data/genomic.fna -o ${date_today}.sequences.zst
    rm tmp.zip
    rm -rf ncbi_dataset
    '''

}

process extract_new_records {
    // Keeps only new records added since previous run
    cpus 1
    penv 'smp'
    conda "python=3.9 zstd seqtk"

    input:
        tuple file(metadata_json), file(sequences_fasta)
    output:
        file '*.tsv'

    shell:
    '''
    date_today=$(date +%Y-%m-%d)

    python !{primer_monitor_path}/lib/parse_ncbi.py <(zstd -d --long=30 < !{metadata_json}) <(zstd -d --long=30 < !{prev_json}) <(zstd -d --long=30 < !{sequences_fasta} | seqtk seq | paste - -) ${date_today}.tsv

    find !{output_path} -maxdepth 1 -mtime +5 -type f -name "*.metadata.zst" -name "*.sequences.zst" -delete
    '''

}

process transform_data {
    cpus 1
    penv 'smp'
    conda "gawk"

    input:
        file ncbi_tsv

    output:
        tuple file('*.metadata'), file('*.fasta')

    shell:
    '''
    date_today=$(date +%Y-%m-%d)

    # Create an empty FASTA in case there are no new seqs
    touch ${date_today}.fasta
    cat !{ncbi_tsv} | gawk -F'\t' -f !{primer_monitor_path}/lib/process_seqs.awk -v cur_date=${date_today}
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

process get_pangolin_version {
    cpus 1

    output:
        env pangolin_version
        env pangolin_data_version
        env use_pending
    shell:
    '''
    attempts=0
    while [ -f "!{flag_path}/pangolin_version_mutex.lock" ]; do
            if [ "$attempts" -gt 10 ]; then
                exit 1
            fi
            sleep 60
            attempts=$((attempts + 1))
    done
    touch !{flag_path}/pangolin_version_mutex.lock;
    touch !{flag_path}/summarize_variants_running.lock;
    use_pending="false"
    if [ -f "!{flag_path}/recall_pangolin_running.lock" ]; then
        use_pending="true"
    fi
    pangolin_version=$(cat !{params.pangolin_version_path})
    pangolin_data_version=$(cat !{params.pangolin_data_version_path})
    rm !{flag_path}/pangolin_version_mutex.lock;
    '''
    }

process pangolin_calls {
    cpus 8
    conda "pangolin=$pangolin_version pangolin-data=$pangolin_data_version"
    input:
        val pangolin_version
        val pangolin_data_version
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
    input:
        file csv
        file complete
        val use_pending
        //the .complete is only here to make sure this happens *after* the main DB load
    output:
        file '*.complete_pangolin'
    shell:
    '''
    field="pangolin_call_id"
    if [ "$use_pending" = "true" ]; then
        field="pending_pangolin_call_id"
    fi
    PGPASSFILE="!{primer_monitor_path}/config/.pgpass" !{primer_monitor_path}/lib/pangolin_calls/update_fasta_records.sh !{csv} $field
    '''
}

process recalculate_database_views {
    cpus 1
    publishDir "${output_path}", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 2
    input:
        file everything
        file everything_pangolin
    output:
        file refresh_complete.txt
    shell:
    '''
    # recalculate all the views at the end to save time
    RAILS_ENV=production ruby /mnt/bioinfo/prg/primer_monitor/upload.rb --skip_data_import && touch refresh_complete.txt
    '''
}

process update_new_calls {
    cpus 1
    penv 'smp'
    input:
        file all_done
    shell:
    '''
    if [ ! -f "!{flag_path}/recall_pangolin_running.lock" ]; then
        PGPASSFILE="!{primer_monitor_path}/config/.pgpass" !{primer_monitor_path}/lib/pangolin_calls/swap_calls.sh; touch done.txt;
    fi
    rm !{flag_path}/summarize_variants_running.lock
    '''
}

workflow {
    download_data()
    extract_new_records(download_data.out)
    transform_data(extract_new_records.out.splitText(file: true, by: 10000).filter{ it.size()>77 })
    align(transform_data.out)
    load_to_db(align.out)
    get_pangolin_version()
    pangolin_calls(get_pangolin_version.out[0], get_pangolin_version.out[1], transform_data.out)
    load_pangolin_data(pangolin_calls.out, load_to_db.out, get_pangolin_version.out[2])
    recalculate_database_views(load_to_db.out.collect(), load_pangolin_data.out.collect())
    update_new_calls(recalculate_database_views.out)
}

workflow.onError {
    println "removing lock files..."
    //get rid of "pipeline running" lock
    running_lock = file('${flag_path}/summarize_variants_running.lock')
    running_lock.delete()
    //get rid of mutex for pangolin version if it exists and is owned by this pipeline
    mutex = file('${flag_path}/pangolin_version_mutex.lock')
    String mutex_content = ""
    if(mutex.exists()) {
        mutex.withReader {
             mutex_content = it.getText()
        }
        if(mutex_content.equals("summarize_variants"))
        {
            mutex.delete()
        }
    }
}
