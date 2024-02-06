nextflow.enable.dsl = 2

assert params.ref != null : "--ref must be specified"
ref = params.ref
ref = file(ref).toAbsolutePath()

assert params.primer_monitor_path != null : "--primer_monitor_path must be specified"
primer_monitor_path = params.primer_monitor_path

assert params.output_path != null : "--output_path must be specified"
output_path = params.output_path

params.temp_dir = '/tmp'
temp_dir = params.temp_dir

assert params.lineage_caller != null : "--lineage_caller must be specified"
lineage_caller = params.lineage_caller

assert params.lineage_caller_script != null : "--lineage_caller_script must be specified"
lineage_caller_script = params.lineage_caller_script

assert params.taxon_id != null : "--taxon_id must be specified"
taxon_id = params.taxon_id

process download_data {
    // Downloads the full dataset
    cpus 16
    conda "ncbi-datasets-cli unzip zstd"
    errorStrategy 'retry'
    maxRetries 2
    publishDir "${output_path}", mode: 'copy', pattern: '*.zst', overwrite: true
    // mode "link" assumes that the output path is on the same disk as the work directory, switch to copy if not
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
    // Keeps only new records added since previous run
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
    <(zstd -d --long=30 < !{metadata_json}) \
    <(zstd -d --long=30 < !{sequences_fasta} | seqtk seq | paste - -) \
    <(psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -c "SELECT genbank_accession FROM fasta_records;" --csv -t) \
    ${date_today}.tsv;

    find !{output_path} -maxdepth 1 -mtime +5 -type f -name "*.metadata.zst" -name "*.sequences.zst" -delete;
    '''

}

process transform_data {
    cpus 1
    conda "gawk"

    input:
        tuple val(index), file(ncbi_tsv)

    output:
        tuple val(index), file('*.metadata'), file('*.fasta')

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
    conda "minimap2=2.17 sed python=3.9 samtools=1.11"
    publishDir "${output_path}", mode: 'copy', pattern: '*.bam', overwrite: true

    input:
        tuple val(index), file(metadata), file(fasta)
    output:
        tuple val(index), file('*.metadata'), file('*.tsv')

    shell:
    '''

        export TMPDIR="!{temp_dir}"

        date_today=$(date +%Y-%m-%d)

        sed -E 's/ /_/g' !{fasta} \
        | minimap2 -r 10000 --score-N=0 -t !{task.cpus} --eqx -x map-ont -a !{ref} /dev/stdin \
        | samtools view -h -F 2304 /dev/stdin \
        | tee >(samtools view -b -o ${date_today}.$(basename $PWD).bam /dev/stdin) \
        | python3 !{primer_monitor_path}/lib/parse_alignments.py > $(basename $PWD).tsv

        mv !{metadata} $(basename $PWD).metadata

    '''
}

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
    -v "caller_name=!{lineage_caller}" <<< "SELECT version_specifiers FROM lineage_callers WHERE name=:'caller_name';" -t --csv);
    '''
    }

process lineage_calls {
    cpus 8
    conda "${version_spec}"
    input:
        val version_spec
        tuple val(index), file(metadata), file(fasta)
    output:
        tuple val(index), file("*.csv")
    shell:
    '''

    export TMPDIR="!{temp_dir}"

    source "!{primer_monitor_path}/.env"
    export BACKEND_INSTALL_PATH

    !{primer_monitor_path}/lib/lineage_calling/caller_wrappers/!{lineage_caller_script}.sh -@ 8 -d !{taxon_id}/current !{fasta}
    '''
}

process load_to_db {
    cpus 1
    publishDir "${output_path}", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 10
    maxForks 1

    conda "'postgresql>=15'"

    input:
        tuple val(index), file(metadata), file(tsv), file(csv)
    output:
        file '*.complete'
    shell:
    '''
    RAILS_ENV=production ruby !{primer_monitor_path}/upload.rb \
            --import_calls \
            --import_seqs \
            --lineage_csv !{csv} \
            --metadata_tsv !{metadata} \
            --variants_tsv !{tsv} \
            --taxon !{taxon_id} \
            --caller !{lineage_caller} \
            && mv !{metadata} !{metadata}.complete
    '''
}

workflow {
    download_data()
    extract_new_records(download_data.out)
    index = 0
    transform_data(extract_new_records.out.splitText(file: true, by: 10000).filter{ it.size()>77 }.map{ [index++, it] })
    align(transform_data.out)
    get_caller_version()
    lineage_calls(get_caller_version.out, transform_data.out)

    //ensure batches of alignments and lineage calls stay in sync for DB load
    seq_recs = align.out.join(lineage_calls.out)

    load_to_db(seq_recs)
}