#!/usr/bin/env bash

if (($# < 1)); then # if bowtie2 index missing
  echo "usage: get_primers.sh <conda env path> <bowtie2 index> [primer_set_ids...]" >&2
  exit 1
fi

if (($# < 2)); then # if primer set IDs missing
  echo "no primer set ids specified, exiting" >&2
  exit 1
fi

conda_env_path="$1"
bt2_index="$2"
shift 2

# you need to export DB_HOST, DB_NAME, and DB_USER before running this

db_csv=$(mktemp)

null_check="AND ref_start IS NULL"

if [ -n "$REALIGN" ]; then
  null_check=""
fi

for id in "$@"; do
  echo "Processing $id..." >&2

  psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -c "SELECT id, sequence FROM oligos WHERE primer_set_id=$id $null_check;" --csv -t | \
  awk 'BEGIN { FS="," }; {print ">" $1 "\n" $2}' | \
  "$MICROMAMBA_BIN_PATH/micromamba" run -p "$conda_env_path" bowtie2 -f --end-to-end --score-min L,-0.6,-1.5 -L 8 -x "$bt2_index" -U - | \
  "$MICROMAMBA_BIN_PATH/micromamba" run -p "$conda_env_path" samtools view -b | \
  "$MICROMAMBA_BIN_PATH/micromamba" run -p "$conda_env_path" bedtools bamtobed -i - | awk '{print $1 "," $2 "," $3 "," $4}' >>"$db_csv"

done

PGPASSFILE="$PGPASSFILE" psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" >&2 <<CMDS
create temporary table tmp_oligo_alignment_positions (ref_name text, ref_start integer, ref_end integer, seq_id integer);
\copy tmp_oligo_alignment_positions from '$db_csv' with (format csv);
delete from oligo_alignment_positions oap where exists (select 1 from tmp_oligo_alignment_positions top where top.seq_id=oap.oligo_id);
insert into oligo_alignment_positions (oligo_id, organism_taxon_id, ref_start, ref_end, created_at, updated_at)
select seq_id, (select id from organism_taxa where organism_taxa.reference_accession=ref_name), ref_start, ref_end, NOW(), NOW()
from tmp_oligo_alignment_positions;
CMDS

rm "$db_csv"
