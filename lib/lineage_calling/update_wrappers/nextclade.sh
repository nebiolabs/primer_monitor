#!/usr/bin/env bash

# pangolin wrapper script

organism_slug="$1"
taxon_id="$2"

while read -d " " -r version; do
  # if this isn't a version, stop
  if ! grep "=" <<< "$version"; then
    break
  fi

  # packages use 2 equals signs here, (non-conda) datasets use only 1
  if grep "==" <<< "$version"; then
    package_name="$(cut -f 1 -d "=" <<< "$version")"

    latest_version=$("$MICROMAMBA_BIN_PATH/micromamba" search -c bioconda "$package_name" | grep -E "Version[[:blank:]]+[0-9]" | awk '{ print $2 }')
    if [ "$latest_version" = "" ]; then
      # skip this entire taxon and email an error
      echo "Error: Got blank version string when trying to update package '$package_name' of caller 'nextclade' \
      for taxon '$taxon_id' of organism '$organism_slug'. Skipping this taxon." | \
      mail -r "$ADMIN_EMAIL" -s "Lineage caller update error (nextclade:$package_name - $organism_slug)" "$NOTIFICATION_EMAILS"
      exit 1;
    fi
    # append version to $new_caller_version
    new_caller_version="$new_caller_version $package_name==$latest_version"
  else
    # this is a version number for a dataset
    dataset_name="$(cut -f 1 -d "=" <<< "$version")"
    # delete the old dataset, then install a new one
    rm -rf "$BACKEND_INSTALL_PATH/datasets/$taxon_id"
    nextclade dataset get --name="$dataset_name" --output-dir="$BACKEND_INSTALL_PATH/datasets/$taxon_id"
    dataset_version="$(python $BACKEND_INSTALL_PATH/lib/lineage_calling/update_wrappers/extract_dataset_version.py \
    "$BACKEND_INSTALL_PATH/datasets/$(basename "$dataset_name")/pathogen.json")"
    success="$?"
    if [ "$success" -ne 0 ]; then
      # skip this entire taxon and email an error
      echo "Error: Failed to update dataset '$dataset_name' of caller 'nextclade' \
      for taxon '$taxon_id' of organism '$organism_slug'. Skipping this taxon." | \
      mail -r "$ADMIN_EMAIL" -s "Dataset update error (nextclade:$dataset_name - $organism_slug)" "$NOTIFICATION_EMAILS"
      exit 1;
    fi
    new_caller_version="$new_caller_version $dataset_name=$dataset_version"
  fi

done
# trailing space is intentional to make sure read gets the last word in the string

# trim leading space
sed -E "s/^ //" <<< "$new_caller_version"