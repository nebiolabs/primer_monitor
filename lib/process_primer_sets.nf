nextflow.enable.dsl = 2

primer_names = params.primer_names
primer_names = file(primer_names).toAbsolutePath()

params.pct_cutoff = 1
pct_cutoff = params.pct_cutoff

params.score_cutoff = 80
score_cutoff = params.score_cutoff

params.primer_monitor_path = '/mnt/bioinfo/prg/primer_monitor'
primer_monitor_path = params.primer_monitor_path

params.organism_dirname = "2697049"
organism_dirname = params.organism_dirname

process recompute_affected_primers {
    cpus 8
    errorStrategy 'retry'
    maxRetries 2
    conda "libiconv psycopg2 bedtools coreutils 'postgresql>=15' gawk"
    input:
        file primer_names_file
    output:
        file 'primers_done.txt'
    shell:
    '''
    # recompute the primer data for igvjs visualization
    !{primer_monitor_path}/lib/visualization/recompute_affected_primers.sh !{primer_monitor_path} !{organism_dirname} \
    !{pct_cutoff} !{score_cutoff} !{task.cpus} !{primer_names_file};

    cp !{primer_names_file} primers_done.txt
    '''
}

process update_db {
    cpus 1
    errorStrategy 'retry'
    maxRetries 2
    conda "'postgresql>=15'"
    input:
        file completed_primers
    shell:
    '''
    while read -f primer_set; do
        psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -v "primer_set=$primer_set" <<< "UPDATE primer_sets SET status='complete' WHERE name=:'primer_set';";
    done < !{completed_primers}
    '''
}


workflow {
    recompute_affected_primers(params.primer_names)
    update_db(recompute_affected_primers.out)
}
