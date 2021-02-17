ref = params.ref 
ref = file(ref).toAbsolutePath()

prev_json_path = Channel.fromPath(params.prev_json)
// prev_json = file(prev_json).toAbsolutePath()

ncov_path = '/mnt/home/mcampbell/src/ncov-ingest'
primer_monitor_path = '/mnt/bioinfo/prg/primer_monitor'

process download_data {
    // Downloads the full dataset, then keeps only new records added since yesterday
    cpus 1
    conda "curl"
    publishDir "output", mode: 'copy', pattern: 'gisaid.sorted_json', overwrite: true

    input:
        file(prev_json) from prev_json_path
    output:
        file('*json') into downloaded_data
        file('gisaid.sorted_json')

    shell:
    '''
    source !{primer_monitor_path}/.env
    mv !{prev_json} !{prev_json}.old
    date_today=$(date +%Y-%m-%d)
    curl -u $GISAID_USER:$GISAID_PASSWORD https://www.epicov.org/epi3/3p/neb/export/provision.json.xz | xz -d -T8 > tmp.json

    sort tmp.json > gisaid.sorted_json
    rm tmp.json

    comm -13 !{prev_json}.old gisaid.sorted_json > ${date_today}.json
    '''

}

process transform_data {
    cpus 1
    conda "regex fsspec pandas typing"

    input:
        file(gisaid_json) from downloaded_data
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
    publishDir "output", mode: 'copy'

    input:
        file(variants) from split_variants.collect()
    output:
        file("combined_variants.tsv") into combined_variants

    shell:
    '''
    cat !{variants} >> combined_variants.tsv
    '''
}
//    echo -e 'Sample\tReference position\tMismatch type\tMismatch' > combined_variants.tsv


process load_to_db {
    cpus 1

    input:
        tuple file(metadata), file(variants) from transformed_metadata.concat(combined_variants).collect()

    shell:
    '''
    ruby !{primer_monitor_path}/upload.rb --metadata_tsv !{metadata} --variants_tsv !{variants}
    '''

}
