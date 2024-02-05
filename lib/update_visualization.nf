nextflow.enable.dsl = 2

params.pct_cutoff =
pct_cutoff = params.pct_cutoff

params.score_cutoff =
score_cutoff = params.score_cutoff

params.primer_monitor_path =
primer_monitor_path = params.primer_monitor_path

params.override_path =
override_path = params.override_path
override_path = file(override_path).toAbsolutePath()

params.temp_dir = '/tmp'
temp_dir = params.temp_dir

params.organism =
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