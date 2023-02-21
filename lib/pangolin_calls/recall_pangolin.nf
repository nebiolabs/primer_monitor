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

process extract_new_records {
    // Get all (deduplicated) records as of yesterday's run
    cpus 1
    penv 'smp'
    conda "python=3.9 zstd"

    output:
    file '*.json'

    shell:
    '''
    date_yesterday=$(date --date="yesterday" +%Y-%m-%d)
    touch known_empty.json
    python3 !{primer_monitor_path}/lib/filter_duplicates.py <(cat known_empty.zst) <(zstd -d --long=30 < !{full_json}) > ${date_yesterday}.json
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
        file '*.complete_pangolin'
    shell:
    '''
    PGPASSFILE="!{primer_monitor_path}/config/.pgpass" !{primer_monitor_path}/lib/pangolin_calls/update_fasta_records.sh !{csv}
    '''
}


workflow {
    extract_new_records()
    transform_data(extract_new_records.out.splitText(file: true, by: 2500))
    pangolin_calls(transform_data.out)
    load_pangolin_data(pangolin_calls.out, load_to_db.out)
}