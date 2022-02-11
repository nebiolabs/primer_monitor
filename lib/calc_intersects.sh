designs="ARTIC_v3.bed ARTIC_v4.bed ARTIC_v4.1.bed ARTIC_v4.1_(alts_only).bed ARTIC_v4.1_(spiked_alts).bed VarSkip_1a.bed VarSkip_2a.bed Varskip_Long_1a.bed Midnight_1200.bed Resende.bed Swift.bed UNZA.bed USydney.bed"
#for f in $designs; do 
#  bedtools intersect -a $f -b voc_variants.sort.bed > $f.overlaps
#done

for var in {AY.,B.,BA.,C.,P.,R.,1000,4000}*.bed ; do 
  for design in $designs; do
    bedtools intersect -wb -a $design -b $var  > ${design}_vs_${var}.overlaps
  done
done

wc -l *.overlaps | head -n -1 | sed -E 's/^ +//' | tr ' ' '\t' > variant_overlap_summary.tsv
for f in *.overlaps; do sed "s/^/$f\t/" $f ; done  > variant_overlaps_combined.tsv

