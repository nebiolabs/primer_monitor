#!/usr/bin/env bash

source "$(dirname "$0")/../../.env";
mkdir -p "$BACKEND_SCRATCH_PATH/status";

"$NEXTFLOW_INSTALL_PATH" -log "$BACKEND_SCRATCH_PATH/log_download/" \
run "$BACKEND_INSTALL_PATH/lib/summarize_variants.nf" \
--ref "$BACKEND_INSTALL_PATH/igvstatic/2697049/ref/NC_045512.2.fasta" \
-w "$BACKEND_SCRATCH_PATH/work_download/" \
--pangolin_version_path "$BACKEND_INSTALL_PATH/pangolin_ver.txt" \
--pangolin_data_version_path "$BACKEND_INSTALL_PATH/pangolin_data_ver.txt" \
--flag_path "$BACKEND_SCRATCH_PATH/status" \
--igvstatic_path "$FRONTEND_HOST" \
--frontend_host "$FRONTEND_HOST" \
--jump_proxy "$JUMP_PROXY" \
-N "$NOTIFICATION_EMAILS";