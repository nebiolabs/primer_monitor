#!/usr/bin/env bash

# This is run by cron at 2am every morning.
# The path to the .env file is passed as the first argument.

dotenv_path=$1

source "$(dirname "$0")/../echo_log.sh"

# shellcheck source=../../.env
source "$dotenv_path";

export PATH="$PATH:$MICROMAMBA_BIN_PATH:$CONDA_BIN_PATH:$QSUB_PATH"
export NXF_CONDA_CACHEDIR="$BACKEND_SCRATCH_PATH/conda_envs"
export NXF_JAVA_HOME

export TMPDIR="${TEMP_DIR:-/tmp}"

echo_log "======= START ======="

while read -r taxon; do
  organism_slug="$(cut -f 1 -d "," <<< "$taxon")"
  ref_accession="$(cut -f 2 -d "," <<< "$taxon")"
  caller_name="$(cut -f 3 -d "," <<< "$taxon")"
  caller_script_name="$(cut -f 4 -d "," <<< "$taxon")"
  taxon_id="$(cut -f 5 -d "," <<< "$taxon")"

  echo_log "processing taxon $taxon_id for $organism_slug, accession $ref_accession, caller $caller_name"

  "$NEXTFLOW_INSTALL_PATH" -quiet -log "$BACKEND_SCRATCH_PATH/log_download-$(date +%F_%T)" \
  run "$BACKEND_INSTALL_PATH/lib/summarize_variants.nf" \
  --ref "$BACKEND_INSTALL_PATH/igvstatic/$organism_slug/ref/$ref_accession.fasta" \
  -w "$BACKEND_SCRATCH_PATH/work_download/" \
  --primer_monitor_path "$BACKEND_INSTALL_PATH" \
  --output_path "$BACKEND_SCRATCH_PATH" \
  --lineage_caller "$caller_name" \
  --lineage_caller_script "$caller_script_name" \
  --taxon_id "$taxon_id" \
  --temp_dir "$TMPDIR" \
  -N "$NOTIFICATION_EMAILS" \
  -c "$BACKEND_INSTALL_PATH/lib/nextflow.config";

  echo_log "processed taxon $taxon_id"

done < <("$PSQL_INSTALL_PATH" -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" \
-c "SELECT o.slug,ot.reference_accession,lc.name,lc.script_name,ot.ncbi_taxon_id \
FROM organisms o INNER JOIN organism_taxa ot ON ot.organism_id=o.id LEFT JOIN lineage_callers lc \
ON ot.lineage_caller_id=lc.id WHERE o.public IS TRUE;" -t --csv);

while read -r organism; do
  organism_slug="$(cut -f 1 -d "," <<< "$organism")"

  echo_log "updating visualization data for $organism_slug"

  "$NEXTFLOW_INSTALL_PATH" -quiet -log "$BACKEND_SCRATCH_PATH/log_visualization-$(date +%F_%T)" \
  run "$BACKEND_INSTALL_PATH/lib/update_visualization.nf" \
  -w "$BACKEND_SCRATCH_PATH/work_visualization/" \
  --primer_monitor_path "$BACKEND_INSTALL_PATH" \
  --pct_cutoff "$PCT_CUTOFF" \
  --score_cutoff "$SCORE_CUTOFF" \
  --override_path "$BACKEND_INSTALL_PATH/igvstatic/$organism_slug/overrides.txt" \
  --organism "$organism_slug" \
  -c "$BACKEND_INSTALL_PATH/lib/nextflow.config";

  echo_log "updated visualization data for $organism_slug"

done < <("$PSQL_INSTALL_PATH" -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" \
-c "SELECT slug FROM organisms WHERE public IS TRUE;" -t --csv);

echo_log "rebuilding materialized views"
PGPASSFILE="$BACKEND_INSTALL_PATH/config/.pgpass" "$PSQL_INSTALL_PATH" -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" \
-c "REFRESH MATERIALIZED VIEW variant_overlaps; \
REFRESH MATERIALIZED VIEW counts; \
REFRESH MATERIALIZED VIEW time_counts; \
REFRESH MATERIALIZED VIEW oligo_variant_overlaps; \
REFRESH MATERIALIZED VIEW identify_primers_for_notifications; \
REFRESH MATERIALIZED VIEW initial_score; \
REFRESH MATERIALIZED VIEW lineage_info;";
echo_log "done"

