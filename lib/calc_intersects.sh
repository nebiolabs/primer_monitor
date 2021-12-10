designs="ARTIC_v3.bed ARTIC_v4.bed VarSkip_1a.bed Varskip1_Long.bed Midnight_1200.bed"
#for f in $designs; do 
#  bedtools intersect -a $f -b voc_variants.sort.bed > $f.overlaps
#done

for var in {AY.,B.,BA.,C.,P.,R.,1000,4000}*.bed ; do 
  for design in $designs; do
    bedtools intersect -a $design -b $var  > ${design}_vs_${var}.overlaps
  done
done

wc -l *.overlaps | head -n -1 | sed -E 's/^ +//' | tr ' ' '\t' > variant_overlap_summary.tsv
for f in *.overlaps; do sed "s/^/$f\t/" $f ; done  > variant_overlaps_combined.tsv

