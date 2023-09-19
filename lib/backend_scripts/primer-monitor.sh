#!/usr/bin/env bash

# This is run by cron at 2am every morning.
# The path to the .env file is passed as the first argument.

dotenv_path=$1

# shellcheck source=../../.env
source "$dotenv_path";

if ! ( set -o noclobber; : > "$BACKEND_SCRATCH_PATH/status/summarize_variants_running.lock" ) &> /dev/null; then
    echo "Another summarize_variants instance is running, aborting..." >&2
    exit 1;
fi

mkdir -p "$BACKEND_SCRATCH_PATH/status";

export PATH="$PATH:$CONDA_BIN_PATH:$QSUB_PATH"
export NXF_CONDA_CACHEDIR="$BACKEND_SCRATCH_PATH/conda_download"
export SGE_QMASTER_PORT SGE_ROOT NXF_JAVA_HOME QSUB_PATH

"$NEXTFLOW_INSTALL_PATH" -log "$BACKEND_SCRATCH_PATH/log_download-$(date +%F_%T)" \
run "$BACKEND_INSTALL_PATH/lib/summarize_variants.nf" \
--ref "$BACKEND_INSTALL_PATH/igvstatic/2697049/ref/NC_045512.2.fasta" \
-w "$BACKEND_SCRATCH_PATH/work_download/" \
--primer_monitor_path "$BACKEND_INSTALL_PATH" \
--output_path "$BACKEND_SCRATCH_PATH" \
--pangolin_version_path "$BACKEND_INSTALL_PATH/pangolin_ver.txt" \
--pangolin_data_version_path "$BACKEND_INSTALL_PATH/pangolin_data_ver.txt" \
--flag_path "$BACKEND_SCRATCH_PATH/status" \
--igvstatic_path "$FRONTEND_HOST" \
--frontend_host "$FRONTEND_HOST" \
--jump_proxy "$JUMP_PROXY" \
--override_path "$BACKEND_INSTALL_PATH/igvstatic/2697049/overrides.txt" \
-N "$NOTIFICATION_EMAILS";

rm "$BACKEND_SCRATCH_PATH/status/summarize_variants_running.lock"