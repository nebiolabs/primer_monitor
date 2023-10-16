#!/usr/bin/env bash

# This is run by cron at 2am every morning.
# The path to the .env file is passed as the first argument.

dotenv_path=$1

# shellcheck source=../../.env
source "$dotenv_path";

# ensure this directory exists
mkdir -p "$BACKEND_SCRATCH_PATH/status";

export PATH="$PATH:$MICROMAMBA_BIN_PATH:$CONDA_BIN_PATH:$QSUB_PATH"
export NXF_CONDA_CACHEDIR="$BACKEND_SCRATCH_PATH/conda_envs"
export NXF_JAVA_HOME

export TMPDIR="${TEMP_DIR:-/tmp}"

"$NEXTFLOW_INSTALL_PATH" -log "$BACKEND_SCRATCH_PATH/log_download-$(date +%F_%T)" \
run "$BACKEND_INSTALL_PATH/lib/summarize_variants.nf" \
--ref "$BACKEND_INSTALL_PATH/igvstatic/2697049/ref/NC_045512.2.fasta" \
-w "$BACKEND_SCRATCH_PATH/work_download/" \
--primer_monitor_path "$BACKEND_INSTALL_PATH" \
--output_path "$BACKEND_SCRATCH_PATH" \
--pangolin_version_path "$BACKEND_INSTALL_PATH/pangolin_ver.txt" \
--pangolin_data_version_path "$BACKEND_INSTALL_PATH/pangolin_data_ver.txt" \
--flag_path "$BACKEND_SCRATCH_PATH/status" \
--pct_cutoff "$PCT_CUTOFF" \
--score_cutoff "$SCORE_CUTOFF" \
--override_path "$BACKEND_INSTALL_PATH/igvstatic/2697049/overrides.txt" \
--organism_dirname "2697049" \
--temp_dir "$TMPDIR" \
-N "$NOTIFICATION_EMAILS";
