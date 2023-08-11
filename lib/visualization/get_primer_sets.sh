#!/usr/bin/env bash

# you need to export DB_HOST, DB_NAME, and DB_USER_RO before running this

set -e

primer_sets_file="$3"

if [ -z "$primer_sets_file" ]; then
  unset primer_sets_file
fi

primer_sets_tmp=$(mktemp)

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -c "SELECT primer_sets.name, organisms.reference_accession, \
oligos.ref_start,oligos.ref_end, oligos.name, '0' AS score, \
COALESCE(oligos.strand, '.') AS strand FROM oligos INNER JOIN primer_sets \
ON oligos.primer_set_id=primer_sets.id INNER JOIN organisms ON organisms.id=primer_sets.organism_id \
WHERE (primer_sets.status='complete' ${primer_sets_file:+"OR primer_sets.status='ready'"}) AND oligos.ref_start IS NOT NULL;" --csv -t | tr "," "\t" > "$primer_sets_tmp"

if [ $# -le 0 ]; then
  # to prevent the rm further down from destroying files, fail if $1 not set
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
  --csv -t | awk -F',' '{ print ">" $1 "\n" $2 }' > "$2/$seq_name.fasta"
done < <(cut -f 1 < "$primer_sets_tmp" | sort | uniq)
printf "}\n"

mkdir -p "$1"

# $1 is checked to be set above, -f to not error if nothing present
rm -f "$1"/*.bed

while read -r seq_rec; do
  seq_name=$("$(dirname "$0")/urlify_name.sh" "$(echo "$seq_rec" | cut -f 1)")
  echo "$seq_rec" | cut -f 2-7 >> "$1/$seq_name.bed"
done < <(sort -k3 "$primer_sets_tmp")

rm "$primer_sets_tmp"



