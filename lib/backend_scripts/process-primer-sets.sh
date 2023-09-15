#!/usr/bin/env bash

if ! ( set -o noclobber; : > "$BACKEND_SCRATCH_PATH/status/recomputing_primers.lock" ) &> /dev/null; then
    echo "Another primer recomputation is running, aborting..." >&2
    exit 1;
fi

source "$BACKEND_INSTALL_PATH/.env";
mkdir -p "$BACKEND_SCRATCH_PATH/status";

export PATH="$PATH:$CONDA_BIN_PATH:$QSUB_PATH"

export NXF_CONDA_CACHEDIR="$BACKEND_SCRATCH_PATH/conda_primers"

new_primers_file="$(mktemp)"

"$PSQL_INSTALL_PATH" -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -c "SELECT name FROM primer_sets WHERE status='created';" -t --csv > "$new_primers_file";

line_count="$(wc -l < "$new_primers_file")"

# only run the nextflow if there are actually primers, and no other conflicting script is running
if [ "$line_count" -gt 0 ] && [ ! -e "$BACKEND_SCRATCH_PATH/status/loading_data.lock" ]; then

  "$NEXTFLOW_INSTALL_PATH" -log "$BACKEND_SCRATCH_PATH/log_primer_sets/" \
  run "$BACKEND_INSTALL_PATH/lib/pangolin_calls/process_primer_sets.nf" \
  -w "$BACKEND_SCRATCH_PATH/work_primer_sets/" \
  --primer_monitor_path "$BACKEND_INSTALL_PATH" \
  --output_path "$BACKEND_SCRATCH_PATH" \
  --pct-cutoff "$PCT_CUTOFF" \
  --score-cutoff "$SCORE_CUTOFF" \
  --primer-names "$new_primers_file" \
  -N "$NOTIFICATION_EMAILS";
fi

rm "$new_primers_file";

rm "$BACKEND_SCRATCH_PATH/status/recomputing_primers.lock"