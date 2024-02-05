#!/usr/bin/env bash

# This is run by cron every 2 weeks at 3am.
# The path to the .env file is passed as the first argument.

dotenv_path=$1

# defines echo_log function
source "$(dirname "$0")/../echo_log.sh"

# shellcheck source=../../.env
source "$dotenv_path";

export PATH="$PATH:$MICROMAMBA_BIN_PATH:$CONDA_BIN_PATH:$QSUB_PATH"
export NXF_CONDA_CACHEDIR="$BACKEND_SCRATCH_PATH/conda_envs"
export NXF_JAVA_HOME

export TMPDIR=${TEMP_DIR:-/tmp}

echo_log "======= START ======="

while read -r taxon; do
  organism_slug="$(cut -f 1 -d "," <<< "$taxon")"
  caller_name="$(cut -f 2 -d "," <<< "$taxon")"
  caller_script_name="$(cut -f 3 -d "," <<< "$taxon")"
  caller_version="$(cut -f 4 -d "," <<< "$taxon")"
  dataset_version="$(cut -f 5 -d "," <<< "$taxon")"
  taxon_id="$(cut -f 6 -d "," <<< "$taxon")"

  new_caller_version=""

  echo_log "Starting $caller_name version update on $taxon_id for $organism_slug"

  while read -d " " -r version; do
    # if this isn't a version, stop
    if ! grep "=" <<< "$version"; then
      break
    fi
    package_name="$(cut -f 1 -d "=" <<< "$version")"
    latest_version=$("$MICROMAMBA_BIN_PATH/micromamba" search -c bioconda "$package_name" | grep -E "Version[[:blank:]]+[0-9]" | awk '{ print $2 }')
    if [ "$latest_version" = "" ]; then
      # skip this entire taxon and email an error
      echo_log "Error: Got blank version string when trying to update package '$package_name' of caller '$caller_name' \
      for taxon '$taxon_id' of organism '$organism_slug'. Skipping this taxon." | \
      mail -r "$ADMIN_EMAIL" -s "Lineage caller update error ($caller_name:$package_name - $organism_slug)" "$NOTIFICATION_EMAILS"
      continue 2;
    fi
    # append version to $new_caller_version
    new_caller_version="$new_caller_version $package_name==$latest_version"
  done < <(echo "$caller_version ")
  # trailing space is intentional to make sure read gets the last word in the string

  # trim leading space
  new_caller_version="$(sed -E "s/^ //" <<< "$new_caller_version")"

  # set pending version
  PGPASSFILE="$BACKEND_INSTALL_PATH/config/.pgpass" "$PSQL_INSTALL_PATH" -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" \
  -v "new_caller_version=$new_caller_version" -v "caller_name=$caller_name" \
  <<< "UPDATE lineage_callers SET pending_version_specifiers=:'new_caller_version' WHERE name=:'caller_name';"

  echo_log "$caller_name version updated to $new_caller_version"

  echo_log "Starting $caller_name dataset update on $taxon_id for $organism_slug"

  "$NEXTFLOW_INSTALL_PATH" -quiet -log "$BACKEND_SCRATCH_PATH/log_update-$(date +%F_%T)" \
  run "$BACKEND_INSTALL_PATH/lib/lineage_calling/dataset_update_nextflows/${caller_script_name}.nf" \
  -w "$BACKEND_SCRATCH_PATH/work_dataset_update/" \
  --primer_monitor_path "$BACKEND_INSTALL_PATH" \
  --organism_slug "$organism_slug" \
  --taxon_id "$taxon_id" \
  --caller_name "$caller_name" \
  --caller_version "$new_caller_version" \
  --dataset_version "$dataset_version" \
  -N "$NOTIFICATION_EMAILS" \
  -c "$BACKEND_INSTALL_PATH/lib/nextflow.config"
  update_success="$?"

  if [ "$update_success" -ne 0 ]; then
    echo_log "dataset update failed, skipping taxon"
    continue # skip this entire taxon if it failed
  fi

  echo_log "$caller_name dataset updated"

  echo_log "processing taxon $taxon_id for $organism_slug, caller $caller_name"

  # run the pipeline
  "$NEXTFLOW_INSTALL_PATH" -quiet -log "$BACKEND_SCRATCH_PATH/log_lineages-$(date +%F_%T)" \
  run "$BACKEND_INSTALL_PATH/lib/lineage_calling/recall_lineages.nf" \
  -w "$BACKEND_SCRATCH_PATH/work_lineages/" \
  --primer_monitor_path "$BACKEND_INSTALL_PATH" \
  --lineage_caller "$caller_name" \
  --lineage_caller_script "$caller_script_name" \
  --taxon_id "$taxon_id" \
  --temp_dir "$TMPDIR" \
  -N "$NOTIFICATION_EMAILS" \
  -c "$BACKEND_INSTALL_PATH/lib/nextflow.config"
  success="$?"

  echo_log "processed taxon $taxon_id"

  if [ "$success" -eq 0 ]; then
      echo_log "processing successful, updating current lineage caller version"
      # update actual version
      PGPASSFILE="$BACKEND_INSTALL_PATH/config/.pgpass" \
      "$PSQL_INSTALL_PATH" -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" -v "caller_name=$caller_name" <<SQL
      UPDATE lineage_callers SET version_specifiers=pending_version_specifiers,
      dataset_versions=pending_dataset_versions WHERE name=:'caller_name';
SQL
      # if there is a dataset in $BACKEND_INSTALL_PATH/datasets
      if [ -d "$BACKEND_INSTALL_PATH/datasets/$taxon_id" ]; then
        # move pending to current (and remove old current in the process)
        flock "$LOCK_PATH/primer_monitor_swap_datasets_$ENV_NAME.lock" \
        bash -c "rm -rf $BACKEND_INSTALL_PATH/datasets/$taxon_id/current && \
        mv $BACKEND_INSTALL_PATH/datasets/$taxon_id/pending $BACKEND_INSTALL_PATH/datasets/$taxon_id/current";
      fi
  fi

done < <("$PSQL_INSTALL_PATH" -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" \
-c "SELECT o.slug,lc.name,lc.script_name,lc.version_specifiers,lc.dataset_versions,ot.ncbi_taxon_id \
FROM organisms o INNER JOIN organism_taxa ot ON ot.organism_id=o.id LEFT JOIN lineage_callers lc \
ON ot.caller_id=lc.id;" -t --csv);

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
  --temp_dir "$TMPDIR" \
  -N "$NOTIFICATION_EMAILS" \
  -c "$BACKEND_INSTALL_PATH/lib/nextflow.config";

  echo_log "updated visualization data for $organism_slug"

done < <("$PSQL_INSTALL_PATH" -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER_RO" \
-c "SELECT slug FROM organisms;" -t --csv);

echo_log "rebuilding materialized views"
PGPASSFILE="$BACKEND_INSTALL_PATH/config/.pgpass" "$PSQL_INSTALL_PATH" -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" \
-c "REFRESH MATERIALIZED VIEW lineage_info;";
echo_log "done"