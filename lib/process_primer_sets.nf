nextflow.enable.dsl = 2

assert params.primer_names != null : "--primer_names must be specified"
primer_names = params.primer_names
primer_names = file(primer_names).toAbsolutePath()

assert params.primer_names != null : "--pct_cutoff must be specified"
pct_cutoff = params.pct_cutoff

assert params.score_cutoff != null : "--score_cutoff must be specified"
score_cutoff = params.score_cutoff

assert params.primer_monitor_path != null : "--primer_monitor_path must be specified"
primer_monitor_path = params.primer_monitor_path

assert params.organism != null : "--organism must be specified"
organism = params.organism

assert params.override_path != null : "--override_path must be specified"
override_path = params.override_path
override_path = file(override_path).toAbsolutePath()

process compute_visualization_data {
    cpus 8
    errorStrategy 'retry'
    maxRetries 2
    conda "libiconv psycopg2 bedtools coreutils 'postgresql>=15' gawk"
    input:
        path primer_names_file
    output:
        path 'primers_done.txt'
    shell:
    '''
    #!/usr/bin/env bash

    set -e

    cp !{primer_names_file} primers_done.txt

    # recompute the primer data for igvjs visualization
    !{primer_monitor_path}/lib/visualization/update_visualization_data.sh -o !{override_path} \
    -p "$(pwd)/!{primer_names_file}" !{primer_monitor_path} !{organism} \
    !{pct_cutoff} !{score_cutoff} !{task.cpus};
    '''
}

process update_db {
    cpus 1
    errorStrategy 'retry'
    maxRetries 2
    conda "'postgresql>=15'"
    input:
        path completed_primers
    shell:
    '''
    #!/usr/bin/env bash

    set -e

    source "!{primer_monitor_path}/.env"

    while read -r primer_set; do
        PGPASSFILE="!{primer_monitor_path}/config/.pgpass" psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" \
        -v "primer_set=$primer_set" <<< "UPDATE primer_sets SET status='complete' WHERE name=:'primer_set';";
    done < !{completed_primers}
    '''
}


workflow {
    compute_visualization_data(primer_names)
    update_db(compute_visualization_data.out)
}
