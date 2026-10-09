#!/usr/bin/env bash

if (($# < 1)); then # if bowtie2 index missing
  echo "usage: update_primers.sh <pixi manifest> <bowtie2 index> [primer_set_ids...]" >&2
  exit 1
fi

if (($# < 2)); then # if primer set IDs missing
  echo "no primer set ids specified, exiting" >&2
  exit 1
fi

pixi_manifest="$1"
bt2_index="$2"
shift 2

log() { echo "$(date '+%F %T') update_primers.sh: $*" >&2; }

# put bowtie2, samtools and bedtools from the locked pixi environment on PATH
# (captured first: eval of a failed command substitution is an empty eval, which succeeds)
pixi_hook=$("${PIXI_BIN_PATH:+$PIXI_BIN_PATH/}pixi" shell-hook --frozen --manifest-path "$pixi_manifest") || {
  log "FAILED setting up the pixi environment from $pixi_manifest"
  exit 1
}
eval "$pixi_hook"

# you need to export DB_HOST, DB_NAME and DB_USER before running this

log "aligning primer sets $* to $bt2_index"

db_csv=$(mktemp)
stages=(psql awk bowtie2 samtools bedtools awk)
failed=0

for id in "$@"; do
  log "processing primer set $id"

  psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" -c "SELECT id, sequence FROM oligos WHERE primer_set_id=$id;" --csv -t | \
  awk 'BEGIN { FS="," }; {print ">" $1 "\n" $2}' | \
  bowtie2 -f --end-to-end --score-min L,-0.6,-1.5 -L 8 -x "$bt2_index" -U - | \
  samtools view -b | \
  bedtools bamtobed -i - | awk '{print $1 "," $2 "," $3 "," $4}' >>"$db_csv"

  # the pipeline's own status is only the last awk's, so check each stage
  statuses=("${PIPESTATUS[@]}")
  for i in "${!statuses[@]}"; do
    if ((statuses[i] != 0)); then
      log "primer set $id: ${stages[i]} exited with ${statuses[i]}"
      failed=1
    fi
  done
done

cat "$db_csv" >&2

# a failed stage may have left partial alignments, so keep the old positions rather than load them
if ((failed)); then
  log "FAILED, database not updated; alignments so far are in $db_csv"
  exit 1
fi

if ! PGPASSFILE="$PGPASSFILE" psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" -v ON_ERROR_STOP=1 --single-transaction >&2 <<CMDS
create temporary table tmp_oligo_alignment_positions (ref_name text, ref_start integer, ref_end integer, seq_id integer);
\copy tmp_oligo_alignment_positions from '$db_csv' with (format csv);
delete from oligo_alignment_positions oap where exists (select 1 from tmp_oligo_alignment_positions top where top.seq_id=oap.oligo_id);
insert into oligo_alignment_positions (oligo_id, organism_taxon_id, ref_start, ref_end, created_at, updated_at)
select seq_id, (select id from organism_taxa where organism_taxa.reference_accession=ref_name), ref_start, ref_end, NOW(), NOW()
from tmp_oligo_alignment_positions;
CMDS
then
  log "FAILED loading alignments into the database; they are in $db_csv"
  exit 1
fi

rm "$db_csv"
log "done"
