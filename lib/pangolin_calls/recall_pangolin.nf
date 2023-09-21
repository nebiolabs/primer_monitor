nextflow.enable.dsl=2

params.primer_monitor_path =
primer_monitor_path = params.primer_monitor_path
params.output_path =
output_path = params.output_path

params.pangolin_version_path =
params.pangolin_data_version_path =

params.temp_dir = '/tmp'
temp_dir = params.temp_dir

process get_new_versions {
    cpus 1
    penv 'smp'
    conda "'bash>=4.1'"

    output:
        env latest_pangolin
        env latest_pangolin_data

    shell:
    '''
    #! /usr/bin/env bash

    source "!{primer_monitor_path}/.env"

    touch "$LOCK_PATH/pangolin_version_mutex.lock"
    # gets a file descriptor for the lock file, opened for writing, and saves its number in $lock_fd
    exec {lock_fd}>"$LOCK_PATH/pangolin_version_mutex.lock"
    flock $lock_fd
    export latest_pangolin=$("$CONDA_BIN_PATH/micromamba" search -c bioconda pangolin | grep -E "Version[[:blank:]]+[0-9]" | awk '{ print $2 }')
    export latest_pangolin_data=$("$CONDA_BIN_PATH/micromamba" search -c bioconda pangolin-data | grep -E "Version[[:blank:]]+[0-9]" | awk '{ print $2 }')

    cp "!{params.pangolin_version_path}" "!{params.pangolin_version_path}.old"
    cp "!{params.pangolin_data_version_path}" "!{params.pangolin_data_version_path}.old"

    printf "$latest_pangolin" > "!{params.pangolin_version_path}"
    printf "$latest_pangolin_data" > "!{params.pangolin_data_version_path}"
    # closes the file descriptor in $lock_fd
    exec {lock_fd}>&-
    '''
}

process download_data {
    // Downloads the full dataset
    cpus 16
    conda "ncbi-datasets-cli unzip zstd"
    errorStrategy 'retry'
    maxRetries 2

    output:
        tuple file('*.metadata.zst'), file('*.sequences.zst')

    shell:
    '''
    date_today=$(date +%Y-%m-%d)
    datasets download virus genome taxon SARS-CoV-2 --complete-only --host human --filename tmp.zip
    unzip tmp.zip
    zstd --long=30 --ultra -22 -T!{task.cpus} ncbi_dataset/data/data_report.jsonl -o ${date_today}.metadata.zst
    zstd --long=30 --ultra -22 -T!{task.cpus} ncbi_dataset/data/genomic.fna -o ${date_today}.sequences.zst
    rm tmp.zip
    rm -rf ncbi_dataset
    '''

}


process extract_new_records {
    // Get all (deduplicated) records currently in the database
    cpus 1
    conda "python=3.9 zstd seqtk 'postgresql>=15'"

    input:
        tuple file(metadata_json), file(sequences_fasta)
    output:
        file '*.tsv'


    shell:
    '''

    TMPDIR="!{temp_dir}"
    export TMPDIR

    source "!{primer_monitor_path}/.env"

    date_today=$(date +%Y-%m-%d)

    python !{primer_monitor_path}/lib/parse_ncbi.py \
    --output-existing \
    <(zstd -d --long=30 < !{metadata_json}) \
    <(psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -c "SELECT genbank_accession FROM fasta_records;" --csv -t) \
    <(zstd -d --long=30 < !{sequences_fasta} | seqtk seq | paste - -) \
    ${date_today}.tsv;
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
    date_today=$(date +%Y-%m-%d)

    # Create an empty FASTA in case there are no seqs
    touch ${date_today}.fasta
    cat !{ncbi_tsv} | gawk -F'\t' -f !{primer_monitor_path}/lib/process_seqs.awk -v cur_date=${date_today}
    '''

}

process pangolin_calls {
    cpus 8
    penv 'smp'
    conda "pangolin==$pangolin_version pangolin-data==$pangolin_data_version"
    input:
        val pangolin_version
        val pangolin_data_version
        tuple file(metadata), file(fasta)
    output:
        file "*.csv"
    shell:
    '''

    TMPDIR="!{temp_dir}"
    export TMPDIR

    !{primer_monitor_path}/lib/pangolin_calls/run_pangolin.sh !{fasta} 8
    '''
}

process load_pangolin_data {
    cpus 1
    penv 'smp'

    conda "'postgresql>=15'"

    input:
        file csv
    output:
        file '*.complete_pangolin'
    shell:
    '''
    RAILS_ENV=production ruby !{primer_monitor_path}/upload.rb \
            --import_calls \
            --pangolin_csv !{csv} \
            --pending true \
            && mv !{csv} !{csv}.$(basename $PWD).complete_pangolin
    '''
}

process update_current_calls {
    cpus 1
    penv 'smp'

    conda "'postgresql>=15'"

    input:
        file everything
    output:
        file 'done.txt'
    shell:
    '''
    PGPASSFILE="!{primer_monitor_path}/config/.pgpass" !{primer_monitor_path}/lib/pangolin_calls/swap_current_calls.sh; touch done.txt;
    '''
}

process update_new_calls {
    cpus 1
    penv 'smp'

    conda "'postgresql>=15'"

    input:
        file all_done
    shell:
    '''
    PGPASSFILE="!{primer_monitor_path}/config/.pgpass" !{primer_monitor_path}/lib/pangolin_calls/swap_new_calls.sh; touch done.txt;
    rm "!{params.pangolin_version_path}.old"
    rm "!{params.pangolin_data_version_path}.old"
    '''
}


workflow {
    get_new_versions()
    download_data()
    extract_new_records(download_data.out)
    transform_data(extract_new_records.out.splitText(file: true, by: 2500).filter{ it.size()>77 })
    pangolin_calls(get_new_versions.out[0], get_new_versions.out[1], transform_data.out)
    load_pangolin_data(pangolin_calls.out)
    update_current_calls(load_pangolin_data.out.collect())
    update_new_calls(update_current_calls.out)
}