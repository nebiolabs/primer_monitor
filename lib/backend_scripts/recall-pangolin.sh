#!/usr/bin/env bash

source "$(dirname "$0")/../../.env";
mkdir -p "$BACKEND_SCRATCH_PATH/status";

"$NEXTFLOW_INSTALL_PATH" -log "$BACKEND_SCRATCH_PATH/log_pangolin/" \
run "$BACKEND_INSTALL_PATH/lib/pangolin_calls/recall_pangolin.nf" \
-w "$BACKEND_SCRATCH_PATH/work_pangolin/" \
--pangolin_version_path "$BACKEND_INSTALL_PATH/pangolin_ver.txt" \
--pangolin_data_version_path "$BACKEND_INSTALL_PATH/pangolin_data_ver.txt" \
--flag_path "$BACKEND_SCRATCH_PATH/status" \
--pct-cutoff "$PCT_CUTOFF" \
--score-cutoff "$SCORE_CUTOFF" \
-N "$NOTIFICATION_EMAILS" \
-c "$BACKEND_INSTALL_PATH/lib/nextflow.config"