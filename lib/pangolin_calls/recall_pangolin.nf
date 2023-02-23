nextflow.enable.dsl=2

ncov_path = '/mnt/home/mcampbell/src/ncov-ingest'
primer_monitor_path = '/mnt/bioinfo/prg/primer_monitor'
output_path = '/mnt/hpc_scratch/primer_monitor'

params.conda_path =
params.pangolin_version_path =
params.pangolin_data_version_path =

process get_new_versions {
    cpus 1
    penv 'smp'

    output:
    env latest_pangolin
    env latest_pangolin_data

    shell:
    '''
    touch !{flag_path}/recall_pangolin_running.txt;
    latest_pangolin=$(!{params.conda_path} search -q -c bioconda pangolin | awk '{ print $2 }' | tail -n 1)

    latest_pangolin_data=$(!{params.conda_path} search -q -c bioconda pangolin-data | awk '{ print $2 }' | tail -n 1)
    '''
}

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
    rm known_empty.json
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
    date_yesterday=$(date --date="yesterday" +%Y-%m-%d)
    !{ncov_path}/bin/transform-gisaid --output-metadata ${date_yesterday}.metadata --output-fasta ${date_yesterday}.fasta --output-additional-info ${date_yesterday}.info ${date_yesterday}.*.json
    '''

}


process pangolin_calls {
    cpus 8
    penv 'smp'
    conda "pangolin=$pangolin_version pangolin-data=$pangolin_data_version"
    input:
        val pangolin_version
        val pangolin_data_version
        tuple file(metadata), file(fasta)
        file pangolin_update_complete
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
    output:
        file '*.complete_pangolin'
    shell:
    '''
    PGPASSFILE="!{primer_monitor_path}/config/.pgpass" !{primer_monitor_path}/lib/pangolin_calls/update_fasta_records.sh !{csv} pending_pangolin_call_id
    '''
}

process update_current_calls {
    cpus 1
    penv 'smp'
    input:
        file everything
    output:
        file 'done.txt'
    shell:
    '''
    PGPASSFILE="!{primer_monitor_path}/config/.pgpass" !{primer_monitor_path}/lib/pangolin_calls/swap_calls.sh NOT; touch done.txt;
    '''
}

process update_pangolin_version_files {
    cpus 1
    penv 'smp'
    input:
        val pangolin_version
        val pangolin_data_version
        file db_update_complete
    output:
        file 'done.txt'
    shell:
    '''
    echo !{pangolin_version} > !{params.pangolin_version_path}
    echo !{pangolin_data_version} > !{params.pangolin_data_version_path}
    touch done.txt;
    '''
}

process update_new_calls {
    cpus 1
    penv 'smp'
    input:
        file all_done
    shell:
    '''
    if [ ! -f "!{flag_path}/summarize_variants_running.txt" ]; then
        PGPASSFILE="!{primer_monitor_path}/config/.pgpass" !{primer_monitor_path}/lib/pangolin_calls/swap_calls.sh; touch done.txt;
    fi
    rm !{flag_path}/recall_pangolin_running.txt;
    '''
}


workflow {
    get_new_versions()
    extract_new_records()
    transform_data(extract_new_records.out.splitText(file: true, by: 2500))
    pangolin_calls(get_new_versions.out[0], get_new_versions.out[1], transform_data.out, update_pangolin.out)
    load_pangolin_data(pangolin_calls.out)
    update_current_calls(load_pangolin_data.out.collect())
    update_pangolin_version_files(get_new_versions.out[0], get_new_versions.out[1], update_current_calls.out)
    update_new_calls(update_pangolin_version_files.out)
}