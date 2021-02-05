ref = params.ref 
ref = file(ref).toAbsolutePath()

ncov_path = '/mnt/home/mcampbell/src/ncov-ingest'
primer_path = '/mnt/home/mcampbell/src/primer_monitor'

process download_data {
    cpus 1
    conda "curl"

    input:
        
    output:
        file('*json') into downloaded_data

    shell:
    '''
    date_today=$(date +%Y-%m-%d)
    curl -u <username>:<password> https://www.epicov.org/epi3/3p/neb/export/provision.json.xz | xz -d -T8 > ${date_today}.json
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
    conda "minimap2=2.17 sed"

    input:
        file(fasta) from transformed_fasta.splitFasta(file: true, by: 10000)
    output:
        file('*.tsv') into split_variants

    shell:
    '''
        sed -E 's/ /_/g' !{fasta} \
        | minimap2 -t !{task.cpus} --eqx -x map-ont -a !{ref} /dev/stdin \
        | python3 !{primer_path}/primer_monitor/lib/parse_alignments.py > $(basename $PWD).tsv

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
    echo -e 'Sample\tReference position\tMismatch type\tMismatch' > combined_variants.tsv
    cat !{variants} >> combined_variants.tsv
    '''
}


process load_to_db {
    cpus 1

    input:
        tuple file(metadata), file(variants) from transformed_metadata.concat(combined_variants).collect()

    shell:
    '''
    ruby !{primer_path}/upload.rb --metadata_tsv !{metadata} --variants_tsv !{variants}
    '''

}
