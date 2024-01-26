#!/usr/bin/env bash

# This is run by cron at 2am every morning.
# The path to the .env file is passed as the first argument.

dotenv_path=$1

# shellcheck source=../../.env
source "$dotenv_path";

export PATH="$PATH:$MICROMAMBA_BIN_PATH:$CONDA_BIN_PATH:$QSUB_PATH"
export NXF_CONDA_CACHEDIR="$BACKEND_SCRATCH_PATH/conda_envs"
export NXF_JAVA_HOME

export TMPDIR="${TEMP_DIR:-/tmp}"

while read -r taxon; do
  organism_slug="$(cut -f 1 -d "," <<< "$taxon")"
  ref_accession="$(cut -f 2 -d "," <<< "$taxon")"
  caller_name="$(cut -f 3 -d "," <<< "$taxon")"
  caller_script_name="$(cut -f 4 -d "," <<< "$taxon")"
  taxon_id="$(cut -f 5 -d "," <<< "$taxon")"

  "$NEXTFLOW_INSTALL_PATH" -log "$BACKEND_SCRATCH_PATH/log_download-$(date +%F_%T)" \
  run "$BACKEND_INSTALL_PATH/lib/summarize_variants.nf" \
  --ref "$BACKEND_INSTALL_PATH/igvstatic/$organism_slug/ref/$ref_accession.fasta" \
  -w "$BACKEND_SCRATCH_PATH/work_download/" \
  --primer_monitor_path "$BACKEND_INSTALL_PATH" \
  --output_path "$BACKEND_SCRATCH_PATH" \
  --lineage_caller "$caller_name" \
  --lineage_caller_script "$caller_script_name" \
  --pct_cutoff "$PCT_CUTOFF" \
  --score_cutoff "$SCORE_CUTOFF" \
  --override_path "$BACKEND_INSTALL_PATH/igvstatic/$organism_slug/overrides.txt" \
  --organism "$organism_slug" \
  --taxon_id "$taxon_id" \
  --temp_dir "$TMPDIR" \
  -N "$NOTIFICATION_EMAILS";

done < <("$PSQL_INSTALL_PATH" -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" \
-c "SELECT o.slug,ot.reference_accession,lc.name,lc.script_name,ot.ncbi_taxon_id \
FROM organisms o INNER JOIN organism_taxa ot ON ot.organism_id=o.id LEFT JOIN lineage_callers lc \
ON ot.caller_id=lc.id;" -t --csv);

RAILS_ENV=production ruby "$BACKEND_INSTALL_PATH/upload.rb" --rebuild_views

