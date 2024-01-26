#!/usr/bin/env bash

# pangolin wrapper script

organism_slug="$1"
taxon_id="$2"

while read -d " " -r version; do
  # if this isn't a version, stop
  if ! grep "=" <<< "$version"; then
    break
  fi
  package_name="$(cut -f 1 -d "=" <<< "$version")"
  latest_version=$("$MICROMAMBA_BIN_PATH/micromamba" search -c bioconda "$package_name" | grep -E "Version[[:blank:]]+[0-9]" | awk '{ print $2 }')
  if [ "$latest_version" = "" ]; then
    # skip this entire taxon and email an error
    echo "Error: Got blank version string when trying to update package '$package_name' of caller 'pangolin' \
    for taxon '$taxon_id' of organism '$organism_slug'. Skipping this taxon." | \
    mail -r "$ADMIN_EMAIL" -s "Lineage caller update error (pangolin:$package_name - $organism_slug)" "$NOTIFICATION_EMAILS"
    exit 1;
  fi
  # append version to $new_caller_version
  new_caller_version="$new_caller_version $package_name==$latest_version"
done
# trailing space is intentional to make sure read gets the last word in the string

# trim leading space
sed -E "s/^ //" <<< "$new_caller_version"