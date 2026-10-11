#!/usr/bin/env bash
#$ -N backfill_taxon
#$ -j y
#$ -pe smp 1

# One-off backfill: SARS-CoV-2 records loaded May 2023 - Feb 2024, before multi-organism support, have
# fasta_records.organism_taxon_id and variant_sites.organism_taxon_id NULL, so every view that joins on taxon
# ignores them. All of them predate the first RSV load (2024-02-26), so they are all SARS-CoV-2.
#
# Usage: qsub -o <log> -v BACKFILL_SCRIPT=<path to this script> backfill_sars_cov2_taxon.sh <.env path> \
#          [other env's download_cron.log...]
#
# Runs as an SGE job, so it can't take the cron flock (that lives in /var/lock on the cron host). Instead it waits
# until the latest nightly run of this environment and of each other log given has logged "done" or "FAILED", and
# it stops before the 2am nightly run, resubmitting itself for 06:00 to carry on. Batched by fasta_record_id with a
# commit per batch and periodic VACUUMs to bound bloat; it only touches rows still NULL, so a rerun resumes.

set -euo pipefail

dotenv_path=$1
shift
other_logs=("$@")

# shellcheck source=../../.env
source "$dotenv_path"
source "$BACKEND_INSTALL_PATH/lib/echo_log.sh" # SGE runs a spooled copy of this script, so not relative to $0

export PGPASSFILE="$BACKEND_INSTALL_PATH/config/.pgpass"
psql_rw() {
  "$PSQL_INSTALL_PATH" -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" -v ON_ERROR_STOP=1 -X -q -At "$@"
}

# true once the last START in a backend log is followed by a line that ends the run: "done" or a fatal FAILED
# (per-taxon and per-organism FAILED lines are logged mid-run, and the run carries on)
run_finished() {
  awk '/======= START =======/ { open = 1 }
       / - done$| - FAILED:| - FAILED to (list|rebuild)/ { open = 0 }
       END { exit open }' "$1"
}

# stop between 01:30 and 06:00 so the nightly 02:00 run (and its view refresh) has the database to itself
in_quiet_window() {
  local hhmm=$((10#$(date +%H%M)))
  (( hhmm >= 130 && hhmm < 600 ))
}

resubmit_and_exit() {
  echo_log "stopping for the nightly run; resubmitting for 06:00"
  qsub -a "$(date +%Y%m%d)0600" -o "$SGE_STDOUT_PATH" -v BACKFILL_SCRIPT="$BACKFILL_SCRIPT" "$BACKFILL_SCRIPT" \
    "$dotenv_path" "${other_logs[@]}"
  exit 0
}

if in_quiet_window; then resubmit_and_exit; fi

for log in "$BACKEND_SCRATCH_PATH/download_cron.log" "${other_logs[@]}"; do
  until run_finished "$log"; do
    echo_log "waiting for the run in $log to finish"
    sleep 300
  done
done

taxon_id=$(psql_rw -c "SELECT id FROM organism_taxa WHERE ncbi_taxon_id = 2697049")
first_rsv_record=$(psql_rw -c "SELECT min(fr.id) FROM fasta_records fr JOIN organism_taxa ot \
  ON ot.id = fr.organism_taxon_id WHERE ot.ncbi_taxon_id <> 2697049")
read -r min_id max_id < <(psql_rw -F ' ' -c "SELECT min(id), max(id) FROM fasta_records WHERE organism_taxon_id IS NULL")

if [[ -z "$min_id" ]]; then
  echo_log "nothing to backfill"
  exit 0
fi
if (( max_id >= first_rsv_record )); then
  echo_log "FAILED: NULL-taxon records reach id $max_id, past the first RSV record $first_rsv_record; refusing"
  exit 1
fi

batch=25000
vacuum_every=20 # batches; about 1M variant_sites rows per batch, so bloat stays near 20M rows between vacuums
echo_log "backfilling taxon $taxon_id for fasta_records $min_id..$max_id in batches of $batch on $(hostname)"

n=0
for ((start = min_id; start <= max_id; start += batch)); do
  if in_quiet_window; then resubmit_and_exit; fi
  end=$((start + batch))
  counts=$(psql_rw -c "SET statement_timeout = 0" -1 \
    -c "WITH v AS (UPDATE variant_sites SET organism_taxon_id = $taxon_id
                   WHERE fasta_record_id >= $start AND fasta_record_id < $end AND organism_taxon_id IS NULL
                   RETURNING 1),
             f AS (UPDATE fasta_records SET organism_taxon_id = $taxon_id
                   WHERE id >= $start AND id < $end AND organism_taxon_id IS NULL
                   RETURNING 1)
        SELECT (SELECT count(*) FROM f) || ' records, ' || (SELECT count(*) FROM v) || ' variants'")
  echo_log "ids $start..$((end - 1)): $counts"
  if (( ++n % vacuum_every == 0 )); then
    echo_log "vacuuming"
    psql_rw -c "SET statement_timeout = 0" -c "VACUUM variant_sites"
  fi
done

echo_log "final vacuum analyze"
psql_rw -c "SET statement_timeout = 0" -c "VACUUM ANALYZE variant_sites" -c "VACUUM ANALYZE fasta_records"

remaining=$(psql_rw -c "SELECT count(*) FROM fasta_records WHERE organism_taxon_id IS NULL")
echo_log "done; $remaining fasta_records still without a taxon"
