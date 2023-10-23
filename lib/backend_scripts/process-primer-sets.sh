#!/usr/bin/env bash

# This is run by cron every minute, and checks the database for any newly-approved primer sets.
# If there are any, it recomputes the igvjs visualization data, otherwise it just exits.
# The path to the .env file is passed as the first argument.

dotenv_path=$1

# shellcheck source=../../.env
source "$dotenv_path";

new_primers_file="$(mktemp -p "$BACKEND_SCRATCH_PATH")"

# select all primer sets with status "processing" for which all oligos are aligned (have a ref_start pos),
# because only aligned primer sets can be processed properly
"$PSQL_INSTALL_PATH" -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" \
-c "SELECT name FROM primer_sets WHERE status='processing' AND NOT EXISTS \
(SELECT 1 FROM oligos WHERE oligos.primer_set_id=primer_sets.id AND oligos.ref_start IS NULL);" -t --csv > "$new_primers_file";

line_count="$(wc -l < "$new_primers_file")"

# only run the nextflow if there are actually primers, and no other conflicting script is running
if [ "$line_count" -gt 0 ]; then

  # ensure this directory exists
  mkdir -p "$BACKEND_SCRATCH_PATH/status";

  export PATH="$PATH:$MICROMAMBA_BIN_PATH:$CONDA_BIN_PATH:$QSUB_PATH"
  export NXF_CONDA_CACHEDIR="$BACKEND_SCRATCH_PATH/conda_envs"
  export NXF_JAVA_HOME

  "$NEXTFLOW_INSTALL_PATH" -log "$BACKEND_SCRATCH_PATH/log_primer_sets-$(date +%F_%T)/" \
  run "$BACKEND_INSTALL_PATH/lib/process_primer_sets.nf" \
  -w "$BACKEND_SCRATCH_PATH/work_primer_sets/" \
  --primer_monitor_path "$BACKEND_INSTALL_PATH" \
  --output_path "$BACKEND_SCRATCH_PATH" \
  --pct_cutoff "$PCT_CUTOFF" \
  --score_cutoff "$SCORE_CUTOFF" \
  --primer_names "$new_primers_file" \
  --override_path "$BACKEND_INSTALL_PATH/igvstatic/2697049/overrides.txt" \
  --organism_dirname "2697049" \
  -N "$NOTIFICATION_EMAILS";
fi

rm "$new_primers_file";
