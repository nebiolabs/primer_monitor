// nextclade wrapper script

nextflow.enable.dsl=2

params.organism_slug =
organism_slug = params.organism_slug

params.primer_monitor_path =
primer_monitor_path = params.primer_monitor_path

params.taxon_id =
taxon_id = params.taxon_id

params.caller_name =
caller_name = params.caller_name

params.caller_version =
caller_version = params.caller_version

params.dataset_version =
dataset_version_str = params.dataset_version

process update_caller_dataset {
    cpus 1

    conda "python=3.9 'postgresql>=15' ${caller_version}"

    shell:
    '''
    #! /usr/bin/env bash

    source "!{primer_monitor_path}/.env"

    set -e

    while read -d " " -r version; do
      # if this is not a version, stop
      if ! grep "=" <<< "$version"; then
        break
      fi

      # this is a version number for a dataset
      dataset_name="$(cut -f 1 -d "=" <<< "$version")"
      # delete the old dataset, then install a new one
      rm -rf "!{primer_monitor_path}/datasets/!{taxon_id}/pending"
      nextclade dataset get --name="$dataset_name" --output-dir="!{primer_monitor_path}/datasets/!{taxon_id}/pending"
      dataset_version="$(python "!{primer_monitor_path}/lib/lineage_calling/dataset_update_nextflows/extract_dataset_version.py" \
      "!{primer_monitor_path}/datasets/!{taxon_id}/pending/pathogen.json")"
      success="$?"
      if [ "$success" -ne 0 ]; then
        # skip this entire taxon and email an error
        echo "Error: Failed to update dataset '$dataset_name' of caller '!{caller_name}' \
        for taxon '!{taxon_id}' of organism '!{organism_slug}'. Skipping this taxon." | \
        mail -r "$ADMIN_EMAIL" -s "Dataset update error (!{caller_name}:$dataset_name - !{organism_slug})" "$NOTIFICATION_EMAILS"
        exit 1;
      fi
      new_dataset_version="$new_dataset_version $dataset_name=$dataset_version"

    done < <(echo "!{dataset_version_str} ")
    # trailing space is intentional to make sure read gets the last word in the string

    # trim leading space
    sed -E "s/^ //" <<< "$new_dataset_version"

    # set pending version
    PGPASSFILE="!{primer_monitor_path}/config/.pgpass" psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" \
    -v "new_dataset_version=$new_dataset_version" -v "caller_name=!{caller_name}" \
    <<< "UPDATE lineage_callers SET pending_dataset_versions=:'new_dataset_version' WHERE name=:'caller_name';"
    '''
    }


workflow {
    update_caller_dataset()
}