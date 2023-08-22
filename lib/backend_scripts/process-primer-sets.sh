#!/usr/bin/env bash

source "$(dirname "$0")/../../.env";

new_primers_file="$(mktemp)"

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" -c "SELECT name FROM primer_sets WHERE status='pending';" -t --csv > "$new_primers_file";

line_count="$(wc -l < "$new_primers_file")"

# only run the nextflow if there are actually primers, and no other conflicting script is running
if [ "$line_count" -gt 0 ] && ! compgen -G "$BACKEND_SCRATCH_PATH/status/recomputing_primers.lock" \
&& ! compgen -G "$BACKEND_SCRATCH_PATH/status/loading_data.lock"; then
  "$NEXTFLOW_INSTALL_PATH" -log "$BACKEND_SCRATCH_PATH/log_primer_sets/" \
  run "$BACKEND_INSTALL_PATH/lib/pangolin_calls/process_primer_sets.nf" \
  -w "$BACKEND_SCRATCH_PATH/work_primer_sets/" \
  --pct-cutoff "$PCT_CUTOFF" \
  --score-cutoff "$SCORE_CUTOFF" \
  --primer-names "$new_primers_file" \
  -N "$NOTIFICATION_EMAILS";
fi

rm "$new_primers_file";