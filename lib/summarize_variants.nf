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

params.pangolin_version_path =
params.pangolin_data_version_path =

params.organism_dirname =
organism_dirname = params.organism_dirname

pangolin_version = file(params.pangolin_version_path).text
pangolin_data_version = file(params.pangolin_data_version_path).text

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
    datasets download virus genome taxon 2697049 --complete-only --host human --filename tmp.zip
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

    date_today=$(date +%Y-%m-%d)

    source "!{primer_monitor_path}/.env"

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

process get_pangolin_version {
    cpus 1

    conda "'bash>=4.1'"

    output:
        env pangolin_version
        env pangolin_data_version
    shell:
    '''
    #! /usr/bin/env bash

    source "!{primer_monitor_path}/.env"

    use_pending="false"
    pangolin_version=$(cat !{params.pangolin_version_path})
    pangolin_data_version=$(cat !{params.pangolin_data_version_path})
    '''
    }

process pangolin_calls {
    cpus 8
    conda "pangolin==$pangolin_version pangolin-data==$pangolin_data_version"
    input:
        val pangolin_version
        val pangolin_data_version
        tuple val(index), file(metadata), file(fasta)
    output:
        tuple val(index), file("*.csv")
    shell:
    '''

    TMPDIR="!{temp_dir}"
    export TMPDIR

    !{primer_monitor_path}/lib/pangolin_calls/run_pangolin.sh !{fasta} 8
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
            --pangolin_csv !{csv} \
            --metadata_tsv !{metadata} \
            --variants_tsv !{tsv} \
            && mv !{metadata} !{metadata}.complete
    '''
}

process recalculate_database_views {
    cpus 1
    publishDir "${output_path}", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 2
    input:
        //these files are to make sure all the load_to_db tasks are done first
        file seq_load_complete
    output:
        file 'refresh_complete.txt';
    shell:
    '''
    # recalculate all the views at the end to save time
    RAILS_ENV=production ruby !{primer_monitor_path}/upload.rb --rebuild_views && touch refresh_complete.txt
    '''
}


process recompute_affected_primers {
    cpus 8
    conda "libiconv psycopg2 bedtools coreutils 'postgresql>=15' gawk bc"
    input:
        file complete
    shell:
    '''
    # recompute the primer data for igvjs visualization
    !{primer_monitor_path}/lib/visualization/recompute_affected_primers.sh -o !{override_path} !{primer_monitor_path} !{organism_dirname} !{pct_cutoff} !{score_cutoff} !{task.cpus}
    '''
}

workflow {
    download_data()
    extract_new_records(download_data.out)
    index = 0
    transform_data(extract_new_records.out.splitText(file: true, by: 10000).filter{ it.size()>77 }.map{ [index++, it] })
    align(transform_data.out)
    get_pangolin_version()
    pangolin_calls(get_pangolin_version.out[0], get_pangolin_version.out[1], transform_data.out)

    //ensure batches of alignments and pangolin calls stay in sync for DB load
    seq_recs = align.out.join(pangolin_calls.out)

    load_to_db(seq_recs)
    recalculate_database_views(load_to_db.out.collect())
    recompute_affected_primers(recalculate_database_views.out)
}