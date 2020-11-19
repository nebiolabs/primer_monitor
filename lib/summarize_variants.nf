params.primer_bed = "../n2_e1_primers.bed"
params.fasta = "../gisaid_hcov-19_2020_06_13_01.fasta"
params.ref = "../NC_045512.2.fasta"

input_fasta = Channel.fromPath(params.fasta)

ref = file(params.ref).toAbsolutePath()
primer_bed = file(params.primer_bed).toAbsolutePath()

process align {
    cpus 16
    publishDir "output"
    conda "minimap2=2.17 sambamba=0.7.1 sed"

    input:
        file(fasta) from input_fasta
    output:
        tuple file('*.bam'),file('*.bai') into combined_bam

    shell:
    '''
        sed -E 's/ /_/g' !{fasta} \
        | minimap2 -t !{task.cpus} -x map-ont -a !{ref} /dev/stdin \
        | sambamba view -S -f bam -l 0 /dev/stdin  \
        | sambamba sort -t !{task.cpus} -o combined.bam /dev/stdin
    '''
}

process split_bam_by_region{
    cpus 16
    publishDir "output"
    conda "sambamba=0.7.1 parallel findutils gawk"

    input: 
        tuple file(bam), file(bai) from combined_bam
    output:
        tuple file('*.bam'), file('*.bai') into bams

    shell:
    '''
        sambamba view combined.bam \
        | /usr/bin/awk 'BEGIN {FS="\\t"}; { split($1,a,"/"); split(a[3],b,"-"); a[2] == "USA" ? file=a[2]"_"substr(b[1],1,2)".sam" : file=a[2]".sam" ; print >> file; close file; }'

        sambamba view -H combined.bam > common_header.txt

        find . -maxdepth 1 -name '*.sam' | parallel -j!{task.cpus} \
        "cat common_header.txt {} | sambamba view -S -f bam -o {/.}.bam /dev/stdin"

        rm *.sam
    '''
}

bams.map { bam,bai -> tuple(bam.simpleName,bam,bai) }.transpose().into{region_bams}

process summarize_variants {
    conda 'samtools gawk sed'

    input: 
        tuple region, file(bam), file(bai) from region_bams
    output: 
        tuple region, file('*.tsv') into region_variant_summaries

    shell:
    '''
        samtools mpileup -f !{ref} -d 1000000 -B -L 1000000 -l !{primer_bed} !{bam} \
        | cut -f 1-5 | gawk -F"\\t" -v OFS="\\t" '{gsub("[N\\\\]\\\\.\\\\$\\\\^]","",$5); print}' \
        | sed -E 's/^/!{region}\\t/' \
        > tmp

        awk -F"\\t" -v OFS="\\t" '{
            if(NR==FNR) {
                primers[$2 ":" $3] = $4
            } else { 
                for (p in primers) { 
                    split(p,loc,":"); 
                    start= loc[1]; end=loc[2]

                    if ($3 >= start && $3 <= end) {
                        name=primers[start":"end]
                        print $1,name,$2,$3,$4,$5,$6
                        break;
                    }
                }
            }
        }' !{primer_bed} tmp > !{region}_variant_summary.tsv 
    '''
}

process aggregate_variants {

    publishDir "output", mode: "copy"

    input: 
        file('*') from region_variant_summaries.flatten().toList()
    output: 
        file('*combined.tsv') into combined_variant_summary
    shell:
    ''' 
        cat *.tsv > variant_summary.combined.tsv
    '''
}
