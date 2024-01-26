#!/usr/bin/env bash

# pangolin wrapper script

organism_slug="$1"
taxon_id="$2"

while read -d " " -r version; do
  # if this isn't a version, stop
  if ! grep "=" <<< "$version"; then
    break
  fi

  # this is a version number for a dataset
  dataset_name="$(cut -f 1 -d "=" <<< "$version")"
  # delete the old dataset, then install a new one
  rm -rf "$BACKEND_INSTALL_PATH/datasets/$taxon_id/pending"
  nextclade dataset get --name="$dataset_name" --output-dir="$BACKEND_INSTALL_PATH/datasets/$taxon_id/pending"
  dataset_version="$(python "$BACKEND_INSTALL_PATH/lib/lineage_calling/update_wrappers/extract_dataset_version.py" \
  "$BACKEND_INSTALL_PATH/datasets/$(basename "$dataset_name")/pending/pathogen.json")"
  success="$?"
  if [ "$success" -ne 0 ]; then
    # skip this entire taxon and email an error
    echo "Error: Failed to update dataset '$dataset_name' of caller 'nextclade' \
    for taxon '$taxon_id' of organism '$organism_slug'. Skipping this taxon." | \
    mail -r "$ADMIN_EMAIL" -s "Dataset update error (nextclade:$dataset_name - $organism_slug)" "$NOTIFICATION_EMAILS"
    exit 1;
  fi
  new_dataset_version="$new_dataset_version $dataset_name=$dataset_version"

done
# trailing space is intentional to make sure read gets the last word in the string

# trim leading space
sed -E "s/^ //" <<< "$new_dataset_version"