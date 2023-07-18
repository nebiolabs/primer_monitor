nextflow.enable.dsl = 2

ref = params.ref
ref = file(ref).toAbsolutePath()
params.prev_json=

params.flag_path='/mnt/hpc_scratch/primer_monitor'

params.pct_cutoff = 1
pct_cutoff = params.pct_cutoff

params.score_cutoff = 100
score_cutoff = params.score_cutoff

prev_json = file(params.prev_json, checkIfExists: true).toAbsolutePath()

params.primer_monitor_path = '/mnt/bioinfo/prg/primer_monitor'
primer_monitor_path = params.primer_monitor_path
params.output_path = '/mnt/hpc_scratch/primer_monitor'
output_path = params.output_path
params.igvstatic_path = '/var/www/igvstatic'
igvstatic_path = params.igvstatic_path

params.pangolin_version_path =
params.pangolin_data_version_path =

params.organism_dirname = "2697049"
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
    // Keeps only new records added since previous run
    cpus 1
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
    publishDir "${output_path}", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 10
    maxForks 1

    conda 'postgresql>=15'

    input:
        tuple file(metadata), file(tsv)
    output:
        file '*.complete'
    shell:
    '''
    touch "!{params.flag_path}/loading_data.lock"
    RAILS_ENV=production ruby !{primer_monitor_path}/upload.rb \
        --import_seqs \
        --metadata_tsv !{metadata} \
        --variants_tsv !{tsv} \
        && mv !{metadata} !{metadata}.complete
    '''
}

process get_pangolin_version {
    cpus 1

    conda "'bash>=4.1'"

    output:
        env pangolin_version
        env pangolin_data_version
        env use_pending
    shell:
    '''
    #! /usr/bin/env bash
    touch "!{params.flag_path}/pangolin_version_mutex.lock"
    # gets a file descriptor for the lock file, opened for writing, and saves its number in $lock_fd
    exec {lock_fd}>"!{params.flag_path}/pangolin_version_mutex.lock"
    flock $lock_fd
    use_pending="false"
    if [ -f "!{params.flag_path}/recall_pangolin_running.lock" ]; then
        use_pending="true"
    fi
    pangolin_version=$(cat !{params.pangolin_version_path})
    pangolin_data_version=$(cat !{params.pangolin_data_version_path})
    # closes the file descriptor in $lock_fd
    exec {lock_fd}>&-
    rm "!{params.flag_path}/pangolin_version_mutex.lock"
    touch "!{params.flag_path}/summarize_variants_running.lock"
    '''
    }

process pangolin_calls {
    cpus 8
    conda "pangolin==$pangolin_version pangolin-data==$pangolin_data_version"
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

    conda "'postgresql>=15'"

    input:
        file csv
        file complete
        val use_pending
        //the .complete is only here to make sure this happens *after* the main DB load
    output:
        file '*.complete_pangolin'
    shell:
    '''
    RAILS_ENV=production ruby !{primer_monitor_path}/upload.rb \
            --import_calls \
            --pangolin_csv !{csv} \
            --pending !{use_pending} \
            && mv !{csv} !{csv}.$(basename $PWD).complete_pangolin
    '''
}

process update_new_calls {
    cpus 1

    conda "'postgresql>=15'"

    input:
        //these files are to make sure all the load_to_db and load_pangolin_data tasks are done first
        file seq_load_complete
        file pangolin_calls_complete
    output:
        file 'done.txt'
    shell:
    '''
    if [ ! -f "!{params.flag_path}/recall_pangolin_running.lock" ]; then
        PGPASSFILE="!{primer_monitor_path}/config/.pgpass" !{primer_monitor_path}/lib/pangolin_calls/swap_new_calls.sh; touch done.txt;
    fi
    rm !{params.flag_path}/summarize_variants_running.lock
    rm "!{params.flag_path}/loading_data.lock"
    '''
}

process recalculate_database_views {
    cpus 1
    publishDir "${output_path}", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 2
    input:
        file calls_updated
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
    publishDir "${igvstatic_path}", mode: 'copy'
    errorStrategy 'retry'
    maxRetries 2
    conda "libiconv psycopg2 bedtools coreutils 'postgresql>=15' gawk"
    input:
        file complete
    output:
        file "${organism_dirname}/lineage_variants"
        file "${organism_dirname}/lineage_sets"
        file "${organism_dirname}/primer_sets"
        file "${organism_dirname}/primer_sets_raw"
        file "${organism_dirname}/misc/lineage_sets.json"
        file "${organism_dirname}/misc/tracks.json"
    shell:
    '''
    # recompute the primer data for igvjs visualization
    set -e # fail on error
    source !{primer_monitor_path}/.env
    export DB_HOST
    export DB_USER
    export DB_NAME
    mkdir -p !{organism_dirname}/misc
    mkdir -p !{organism_dirname}/lineage_sets
    mkdir -p !{organism_dirname}/primer_sets_raw
    mkdir -p !{organism_dirname}/lineage_sets

    !{primer_monitor_path}/lib/visualization/get_lineage_data.sh > lineages.csv

    !{primer_monitor_path}/lib/visualization/get_primer_sets.sh !{organism_dirname}/primer_sets_raw > !{organism_dirname}/misc/tracks.json

    psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" -c "SELECT COALESCE(date_collected, date_submitted), COUNT(*) \
    FROM fasta_records GROUP BY COALESCE(date_collected, date_submitted);" --csv -t > seq_counts.csv

    curl https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json \
    | python !{primer_monitor_path}/lib/visualization/get_lineages_to_show.py A,B lineages.csv seq_counts.csv !{organism_dirname}/lineage_sets

    cat <(printf "{") <(ls !{organism_dirname}/lineage_sets \
    | sed -E 's/^(.*)\\.txt$/"\\1": "\\1.*",/') <(echo '"all": "All"}') > !{organism_dirname}/misc/lineage_sets.json

    ls !{organism_dirname}/primer_sets_raw > primer_sets_data.txt

    cat <(ls !{organism_dirname}/lineage_sets) <(echo "all.txt") \
    | xargs !{primer_monitor_path}/lib/visualization/process_primer_sets_with_lineages.sh - "./!{organism_dirname}" !{pct_cutoff} !{score_cutoff} \
    primer_sets_data.txt "./!{organism_dirname}" !{task.cpus}

    rm primer_sets_data.txt
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
    update_new_calls(load_to_db.out.collect(), load_pangolin_data.out.collect())
    recalculate_database_views(update_new_calls.out)
    recompute_affected_primers(recalculate_database_views.out)
}

workflow.onError {
    println "removing lock files..."
    //if it started loading the new seqs, it's not possible to automatically recover
    if(!(file("${params.flag_path}/loading_data.lock").exists()))
    {
        //get rid of "pipeline running" lock
        running_lock = file('${params.flag_path}/summarize_variants_running.lock')
        running_lock.delete()
    }
}
