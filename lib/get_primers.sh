#!/usr/bin/env bash

if (($# < 2)); then # if N variants file provided
    echo "usage: get_primers.sh <bt2_index_dir> <primer_set_id>" >&2;
    exit 1;
fi
# you need to export DB_HOST, DB_NAME, and DB_USER before running this

primers_file=$(mktemp);
bam_file=$(mktemp);
csv_file=$(mktemp);

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" -c "SELECT id, sequence FROM oligos WHERE primer_set_id=$2;" --csv -t | tr "," "\t" | awk '{print ">" $1 "\n" $2}' > "$primers_file";

bowtie2 -f --end-to-end -x "$1" -U "$primers_file" | samtools view -b > "$bam_file";

bedtools bamtobed -i "$bam_file" | awk '{print $4 "," $2 "," $3}' > "$csv_file";

while read -r primer; do
    seq_id=$(echo "$primer" | cut -f 1 -d ",")
    ref_start=$(echo "$primer" | cut -f 2 -d ",")
    ref_end=$(echo "$primer" | cut -f 3 -d ",")
    echo "UPDATE oligos SET (ref_start, ref_end) = ($ref_start, $ref_end) WHERE id=$seq_id;" #echoing only right now for debug
done < "$csv_file"

rm "$primers_file";
rm "$bam_file";
rm "$csv_file";