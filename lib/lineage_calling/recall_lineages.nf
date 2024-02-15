nextflow.enable.dsl=2

assert params.primer_monitor_path != null : "--primer_monitor_path must be specified"
primer_monitor_path = params.primer_monitor_path

assert params.taxon_id != null : "--taxon_id must be specified"
taxon_id = params.taxon_id

assert params.lineage_caller != null : "--lineage_caller must be specified"
lineage_caller = params.lineage_caller

assert params.lineage_caller_script != null : "--lineage_caller_script must be specified"
lineage_caller_script = params.lineage_caller_script

params.temp_dir = '/tmp'
temp_dir = params.temp_dir

process get_caller_version {
    cpus 1

    conda "'bash>=4.1' 'postgresql>=15'"

    output:
        env version_spec
    shell:
    '''
    #! /usr/bin/env bash

    source "!{primer_monitor_path}/.env"

    version_spec=$(PGPASSFILE="!{primer_monitor_path}/config/.pgpass" psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" \
    -v "caller_name=!{lineage_caller}" <<< "SELECT pending_version_specifiers FROM lineage_callers WHERE name=:'caller_name';" -t --csv);
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

    export TMPDIR="!{temp_dir}"

    date_today=$(date +%Y-%m-%d)
    datasets download virus genome taxon !{taxon_id} --complete-only --host human --filename tmp.zip
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

    export TMPDIR="!{temp_dir}"

    source "!{primer_monitor_path}/.env"

    date_today=$(date +%Y-%m-%d)

    python !{primer_monitor_path}/lib/parse_ncbi.py \
    --output-existing \
    <(zstd -d --long=30 < !{metadata_json}) \
    <(zstd -d --long=30 < !{sequences_fasta} | seqtk seq | paste - -) \
    <(psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -c "SELECT genbank_accession FROM fasta_records;" --csv -t) \
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

process lineage_calls {
    cpus 8
    conda "${version_spec}"
    input:
        val version_spec
        tuple file(metadata), file(fasta)
    output:
        file "*.csv"
    shell:
    '''

    export TMPDIR="!{temp_dir}"

    source "!{primer_monitor_path}/.env"
    export BACKEND_INSTALL_PATH

    !{primer_monitor_path}/lib/lineage_calling/caller_wrappers/!{lineage_caller_script}.sh -@ 8 -d !{taxon_id}/pending !{fasta}
    '''
}

process load_lineage_data {
    cpus 1
    errorStrategy 'retry'
    maxRetries 10
    // don't use all the connections
    maxForks 10

    conda "'postgresql>=15'"

    input:
        file csv
    output:
        file '*.complete_lineages'
    shell:
    '''
    RAILS_ENV=production ruby !{primer_monitor_path}/upload.rb \
            --import_calls \
            --lineage_csv !{csv} \
            --pending \
            --taxon !{taxon_id} \
            --caller !{lineage_caller} \
            && mv !{csv} !{csv}.$(basename $PWD).complete_lineages
    '''
}

process update_calls {
    cpus 1
    penv 'smp'

    conda "'postgresql>=15'"

    input:
        file everything
    output:
        file 'done.txt'
    shell:
    '''
    PGPASSFILE="!{primer_monitor_path}/config/.pgpass" !{primer_monitor_path}/lib/lineage_calling/swap_calls.sh !{taxon_id}; touch done.txt;
    '''
}

process cleanup_old_calls {
    cpus 1
    penv 'smp'

    conda "'postgresql>=15'"

    input:
        file update_done
    shell:
    '''
    PGPASSFILE="!{primer_monitor_path}/config/.pgpass" !{primer_monitor_path}/lib/lineage_calling/cleanup_old_calls.sh !{taxon_id};
    '''
}


workflow {
    get_caller_version()
    download_data()
    extract_new_records(download_data.out)
    transform_data(extract_new_records.out.splitText(file: true, by: 2500).filter{ it.size()>77 })
    lineage_calls(get_caller_version.out, transform_data.out)
    load_lineage_data(lineage_calls.out)
    update_calls(load_lineage_data.out.collect())
    cleanup_old_calls(update_calls.out)
}