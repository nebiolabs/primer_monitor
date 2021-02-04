primer_bed = params.primer_bed //= "../n2_e1_primers.bed"
fasta = params.fasta //= "../gisaid_hcov-19_2020_06_13_01.fasta"
ref = params.ref //= "../NC_045512.2.fasta"

input_fasta = Channel.fromPath(fasta).splitFasta(file: true, by: 10000)

ref = file(ref).toAbsolutePath()
primer_bed = file(primer_bed).toAbsolutePath()

process align {
    cpus 16
    conda "minimap2=2.17 sed"

    input:
        file(fasta) from input_fasta
    output:
        file('*.tsv') into split_variants

    shell:
    '''
        sed -E 's/ /_/g' !{fasta} \
        | minimap2 -t !{task.cpus} --eqx -x map-ont -a !{ref} /dev/stdin \
        | python3 /mnt/home/mcampbell/src/primer_monitor/lib/parse_alignments.py > $(basename $PWD).tsv

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
