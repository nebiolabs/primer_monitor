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
    latest_pangolin=$(!{params.conda_path} search -q -c bioconda pangolin | awk '{ print $2 }' | tail -n 1)
    latest_pangolin_data=$(!{params.conda_path} search -q -c bioconda pangolin-data | awk '{ print $2 }' | tail -n 1)

    printf "$latest_pangolin" > !{params.pangolin_version_path}
    printf "$latest_pangolin_data" > !{params.pangolin_data_version_path}

    touch !{flag_path}/recall_pangolin_running.txt;
    '''
}


process extract_new_records {
    // Get all (deduplicated) records as of yesterday's run
    cpus 1
    conda "python=3.9 zstd seqtk"

    output:
        file '*.tsv'


    shell:
    '''
    date_yesterday=$(date --date="yesterday" +%Y-%m-%d)
    touch known_empty.json
    python !{primer_monitor_path}/lib/parse_ncbi.py <(zstd -d --long=30 < !{output_path}/${date_yesterday}.metadata.zst) known_empty.json <(zstd -d --long=30 < !{output_path}/${date_yesterday}.sequences.zst | seqtk seq | paste - -) ${date_yesterday}.tsv
    rm known_empty.zst
    '''

}


process transform_data {
    cpus 1
    conda "gawk"

    input:
        file ncbi_tsv
    output:
        tuple file('*.metadata'), file('*.fasta')


    shell:
    '''
    date_yesterday=$(date --date="yesterday" +%Y-%m-%d)

    # Create an empty FASTA in case there are no seqs
    touch ${date_yesterday}.fasta
    cat !{ncbi_tsv} | gawk -F'\t' -f !{primer_monitor_path}/lib/process_seqs.awk -v cur_date=${date_yesterday}
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
    transform_data(extract_new_records.out.splitText(file: true, by: 2500).filter{ it.size()>77 })
    pangolin_calls(get_new_versions.out[0], get_new_versions.out[1], transform_data.out)
    load_pangolin_data(pangolin_calls.out)
    update_current_calls(load_pangolin_data.out.collect())
    update_new_calls(update_current_calls.out)
}