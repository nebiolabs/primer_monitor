ref = params.ref 
ref = file(ref).toAbsolutePath()

// prev_json_path = Channel.fromPath(params.prev_json)
prev_json = file(params.prev_json).toAbsolutePath()

ncov_path = '/mnt/home/mcampbell/src/ncov-ingest'
primer_monitor_path = '/mnt/bioinfo/prg/primer_monitor'

process download_data {
    // Downloads the full dataset
    cpus 1
    conda "curl"

    output:
        file('gisaid.full_json') into downloaded_data

    shell:
    '''
    source !{primer_monitor_path}/.env
    curl -u $GISAID_USER:$GISAID_PASSWORD https://www.epicov.org/epi3/3p/neb/export/provision.json.xz | xz -d -T8 > gisaid.full_json
    '''

}

process filter_data {
    // Keeps only new records added since previous run
    cpus 1
    conda "python=3.9"
    publishDir "/mnt/flash_scratch/seq-shepherd/", mode: 'copy', pattern: 'gisaid.full_json', overwrite: true

    input:
        // file(prev_json) from prev_json_path
        file(full_json) from downloaded_data
    output:
        file('*.json') into filtered_data
        file('gisaid.full_json')


    shell:
    '''
    date_today=$(date +%Y-%m-%d)

    python3 !{primer_monitor_path}/lib/filter_duplicates.py !{prev_json} !{full_json} > ${date_today}.json
    '''

}

process transform_data {
    cpus 1
    conda "regex fsspec pandas typing"

    input:
        file(gisaid_json) from filtered_data
    output:
        file('*.metadata') into transformed_metadata
        file('*.fasta') into transformed_fasta


    shell:
    '''
    date_today=$(date +%Y-%m-%d)

    !{ncov_path}/bin/transform-gisaid --output-metadata ${date_today}.metadata --output-fasta ${date_today}.fasta --output-additional-info ${date_today}.info ${date_today}.json 
    '''

}

process align {
    cpus 16
    conda "minimap2=2.17 sed python=3.9"

    input:
        file(fasta) from transformed_fasta.splitFasta(file: true, by: 10000)
    output:
        file('*.tsv') into split_variants

    shell:
    '''
        sed -E 's/ /_/g' !{fasta} \
        | minimap2 -t !{task.cpus} --eqx -x map-ont -a !{ref} /dev/stdin \
        | python3 !{primer_monitor_path}/lib/parse_alignments.py > $(basename $PWD).tsv

    '''
}

process combine_variants {
    cpus 1
    publishDir "/mnt/flash_scratch/seq-shepherd/", mode: 'copy'

    input:
        file(variants) from split_variants.collect()
    output:
        file("combined_variants.tsv") into combined_variants

    shell:
    '''
    cat !{variants} >> combined_variants.tsv
    '''
}


process load_to_db {
    cpus 1

    input:
        tuple file(metadata), file(variants) from transformed_metadata.concat(combined_variants).collect()

    shell:
    '''
    ruby !{primer_monitor_path}/upload.rb --metadata_tsv !{metadata} --variants_tsv !{variants}
    '''

}
