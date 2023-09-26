#!/usr/bin/env bash

# This is run by cron every 2 weeks at 3am.
# The path to the .env file is passed as the first argument.

dotenv_path=$1

# shellcheck source=../../.env
source "$dotenv_path";

# ensure this directory exists
mkdir -p "$BACKEND_SCRATCH_PATH/status";

export PATH="$PATH:$MICROMAMBA_BIN_PATH:$CONDA_BIN_PATH:$QSUB_PATH"
export NXF_CONDA_CACHEDIR="$BACKEND_SCRATCH_PATH/conda_envs"
export NXF_JAVA_HOME

export TMPDIR=${TEMP_DIR:-/tmp}

# roll back pangolin update and exit
rollback_pangolin() {
  mv "$BACKEND_INSTALL_PATH/pangolin_ver.txt.old" "$BACKEND_INSTALL_PATH/pangolin_ver.txt";
  mv "$BACKEND_INSTALL_PATH/pangolin_data_ver.txt.old" "$BACKEND_INSTALL_PATH/pangolin_data_ver.txt";
  exit 1;
}

latest_pangolin=$("$CONDA_BIN_PATH/micromamba" search -c bioconda pangolin | grep -E "Version[[:blank:]]+[0-9]" | awk '{ print $2 }')
latest_pangolin_data=$("$CONDA_BIN_PATH/micromamba" search -c bioconda pangolin-data | grep -E "Version[[:blank:]]+[0-9]" | awk '{ print $2 }')

if [ "$latest_pangolin" = "" ] || [ "$latest_pangolin_data" = "" ]; then
  # if either version is blank, fail
  exit 1;
fi

# back up the old pangolin versions before updating
cp "$BACKEND_INSTALL_PATH/pangolin_ver.txt" "$BACKEND_INSTALL_PATH/pangolin_ver.txt.old"
cp "$BACKEND_INSTALL_PATH/pangolin_data_ver.txt" "$BACKEND_INSTALL_PATH/pangolin_data_ver.txt.old"

echo "Starting Pangolin version update"

# roll back the pangolin version update before exiting if something fails or it receives sigterm
trap rollback_pangolin ERR SIGTERM;

# update the pangolin version
printf "%s" "$latest_pangolin" > "$BACKEND_INSTALL_PATH/pangolin_ver.txt"
printf "%s" "$latest_pangolin_data" > "$BACKEND_INSTALL_PATH/pangolin_data_ver.txt"

echo "Pangolin version updated to $latest_pangolin (data $latest_pangolin_data)"

# run the pipeline
"$NEXTFLOW_INSTALL_PATH" -log "$BACKEND_SCRATCH_PATH/log_pangolin-$(date +%F_%T)" \
run "$BACKEND_INSTALL_PATH/lib/pangolin_calls/recall_pangolin.nf" \
-w "$BACKEND_SCRATCH_PATH/work_pangolin/" \
--primer_monitor_path "$BACKEND_INSTALL_PATH" \
--output_path "$BACKEND_SCRATCH_PATH" \
--pangolin_version_path "$BACKEND_INSTALL_PATH/pangolin_ver.txt" \
--pangolin_data_version_path "$BACKEND_INSTALL_PATH/pangolin_data_ver.txt" \
--flag_path "$BACKEND_SCRATCH_PATH/status" \
--pct-cutoff "$PCT_CUTOFF" \
--score-cutoff "$SCORE_CUTOFF" \
--temp_dir "$TMPDIR" \
-N "$NOTIFICATION_EMAILS" \
-c "$BACKEND_INSTALL_PATH/lib/nextflow.config"

# unset the trap since the update was successful
trap - ERR SIGTERM;

# since the pipeline succeeded, remove the .old files
rm "$BACKEND_INSTALL_PATH/pangolin_ver.txt.old"
rm "$BACKEND_INSTALL_PATH/pangolin_data_ver.txt.old"