nextflow.enable.dsl = 2

assert params.primer_names != null : "--pct_cutoff must be specified"
pct_cutoff = params.pct_cutoff

assert params.score_cutoff != null : "--score_cutoff must be specified"
score_cutoff = params.score_cutoff

assert params.primer_monitor_path != null : "--primer_monitor_path must be specified"
primer_monitor_path = params.primer_monitor_path

assert params.override_path != null : "--override_path must be specified"
override_path = params.override_path
override_path = file(override_path).toAbsolutePath()

params.temp_dir = '/tmp'
temp_dir = params.temp_dir

assert params.organism != null : "--organism must be specified"
organism = params.organism

process update_visualization_data {
    cpus 8
    conda "libiconv psycopg2 bedtools coreutils 'postgresql>=15' gawk bc"
    shell:
    '''
    # recompute the primer data for igvjs visualization
    !{primer_monitor_path}/lib/visualization/update_visualization_data.sh -o !{override_path} !{primer_monitor_path} \
    !{organism} !{pct_cutoff} !{score_cutoff} !{task.cpus}
    '''
}

workflow {
    update_visualization_data()
}