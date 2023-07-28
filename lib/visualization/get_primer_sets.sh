#!/usr/bin/env bash

# you need to export DB_HOST, DB_NAME, and DB_USER_RO before running this

set -e

primer_sets_tmp=$(mktemp)

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -c "SELECT primer_sets.name, organisms.reference_accession, \
oligos.ref_start,oligos.ref_end, regexp_replace(oligos.name, '\s', '_', 'g'), '0' AS score, \
COALESCE(oligos.strand, '.') AS strand FROM oligos INNER JOIN primer_sets \
ON oligos.primer_set_id=primer_sets.id INNER JOIN organisms ON organisms.id=primer_sets.organism_id \
WHERE primer_sets.status='complete' AND oligos.ref_start IS NOT NULL;" --csv -t | tr "," "\t" > "$primer_sets_tmp"

if [ $# -le 0 ]; then
  # to prevent the rm further down from destroying files, fail if $1 not set
  exit 1;
fi


urlify_name ()
{
  echo "$1" | iconv -t ASCII//TRANSLIT | sed -E "s/[ \/]/_/g; s/['\"]//g"
}

printf "{"
first=1
while read -r seq_name_raw; do
  seq_name=$(urlify_name "$seq_name_raw")
  if [ $first -eq 0 ]; then
    printf ',\n'
  fi
  printf "\"%s\": \"%s\"" "$seq_name" "$seq_name_raw"
  first=0
done < <(cut -f 1 < "$primer_sets_tmp" | sort | uniq)
printf "}\n"

mkdir -p "$1"

# $1 is checked to be set above, -f to not error if nothing present
rm -f "$1"/*.bed

while read -r seq_rec; do
  seq_name=$(urlify_name "$(echo "$seq_rec" | cut -f 1)")
  echo "$seq_rec" | cut -f 2-7 >> "$1/$seq_name.bed"
done < <(sort -k3 "$primer_sets_tmp")

rm "$primer_sets_tmp"



