nextflow.enable.dsl = 2

ref = params.ref
ref = file(ref).toAbsolutePath()

params.pct_cutoff =
pct_cutoff = params.pct_cutoff

params.score_cutoff =
score_cutoff = params.score_cutoff

params.primer_monitor_path =
primer_monitor_path = params.primer_monitor_path

params.output_path =
output_path = params.output_path

params.override_path =
override_path = params.override_path
override_path = file(override_path).toAbsolutePath()

params.temp_dir = '/tmp'
temp_dir = params.temp_dir

params.lineage_caller =
lineage_caller = params.lineage_caller

params.lineage_caller_script =
lineage_caller_script = params.lineage_caller_script

params.taxon_id =
taxon_id = params.taxon_id

params.organism =
organism = params.organism

process download_data {
    // Downloads the full dataset
    cpus 16
    conda "ncbi-datasets-cli unzip zstd"
    errorStrategy 'retry'
    maxRetries 2
    publishDir "${output_path}", mode: 'link', pattern: '*.zst', overwrite: true
    // mode "link" assumes that the output path is on the same disk as the work directory, switch to copy if not
    output:
        tuple file('*.metadata.zst'), file('*.sequences.zst')

    shell:
    '''

    TMPDIR="!{temp_dir}"
    export TMPDIR

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

    TMPDIR="!{temp_dir}"
    export TMPDIR

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

        TMPDIR="!{temp_dir}"
        export TMPDIR

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

    export PGPASSFILE="!{primer_monitor_path}/config/.pgpass"

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

    TMPDIR="!{temp_dir}"
    export TMPDIR

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

process update_visualization_data {
    cpus 8
    conda "libiconv psycopg2 bedtools coreutils 'postgresql>=15' gawk bc"
    input:
        file complete
    shell:
    '''
    # recompute the primer data for igvjs visualization
    !{primer_monitor_path}/lib/visualization/update_visualization_data.sh -o !{override_path} !{primer_monitor_path} \
    !{organism} !{pct_cutoff} !{score_cutoff} !{task.cpus} !{taxon_id}
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
    update_visualization_data(load_to_db.out.collect())
}