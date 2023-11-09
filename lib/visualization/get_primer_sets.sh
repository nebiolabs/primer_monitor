#!/usr/bin/env bash

# gets primer sets from database

# you need to export DB_HOST, DB_NAME, and DB_USER_RO before running this

set -e

# Usage: get_primer_sets.sh <primer_bed_dir> <primer_fasta_dir> [primer_sets_file]

# primer_bed_dir and primer_fasta_dir are where primer set BED and FASTA files are output
# primer_sets_file is for processing only specific primer sets

primer_bed_dir="$1"
primer_fasta_dir="$2"
primer_sets_file="$3"

if [ -z "$primer_sets_file" ]; then
  unset primer_sets_file
fi

primer_sets_tmp=$(mktemp)

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -c "SELECT primer_sets.name, organism_taxa.reference_accession, \
oligos.ref_start,oligos.ref_end, regexp_replace(oligos.name, '\s', '_', 'g'), '0' AS score, \
COALESCE(oligos.strand, '.') AS strand FROM oligos INNER JOIN primer_sets \
ON oligos.primer_set_id=primer_sets.id \
INNER JOIN organism_taxa ON organism_taxa.id=oligos.organism_taxon_id \
WHERE (primer_sets.status='complete' ${primer_sets_file:+"OR primer_sets.status='processing'"}) AND oligos.ref_start IS NOT NULL;" --csv -t | tr "," "\t" > "$primer_sets_tmp"

if [ $# -le 0 ]; then
  # to prevent the rm further down from destroying files, fail if $primer_bed_dir not set
  exit 1;
fi

printf "{"
first=1
while read -r seq_name_raw; do
  seq_name=$("$(dirname "$0")/urlify_name.sh" "$seq_name_raw")
  if [ $first -eq 0 ]; then
    printf ',\n'
  fi
  printf "\"%s\": \"%s\"" "$seq_name" "$seq_name_raw"
  first=0

  psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -v "primer_set_name=$seq_name_raw" \
  <<< "SELECT regexp_replace(name, '\s', '_', 'g'),sequence FROM oligos WHERE primer_set_id=(SELECT id FROM primer_sets WHERE name = :'primer_set_name');" \
  --csv -t | awk -F',' '{ print ">" $1 "\n" $2 }' > "$primer_fasta_dir/$seq_name.fasta"
done < <(cut -f 1 < "$primer_sets_tmp" | sort | uniq)
printf "}\n"

mkdir -p "$primer_bed_dir"

# $primer_bed_dir is checked to be set above, -f to not error if nothing present
rm -f "$primer_bed_dir"/*.bed

while read -r seq_rec; do
  seq_name=$("$(dirname "$0")/urlify_name.sh" "$(echo "$seq_rec" | cut -f 1)")
  echo "$seq_rec" | cut -f 2-7 >> "$primer_bed_dir/$seq_name.bed"
done < <(sort -k3 "$primer_sets_tmp")

rm "$primer_sets_tmp"



