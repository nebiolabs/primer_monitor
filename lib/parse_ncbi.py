"""
Filters duplicate sequences out of an NCBI dataset given a list of sequences in the database.
Usage: ./parse_ncbi.py <metadata> <old accessions> <new sequences> <output filename>
"""
import sys
import json


def get_if_exists(cur_dict, *keys):
    """
    returns the value from within the dict, or "" if a key is missing
    """
    try:
        for key in keys:
            cur_dict = cur_dict[key]
    except KeyError:
        return ""
    return cur_dict


# maps accession numbers to data
output_lines = {}

# accession numbers seen in prev
prev_accessions = set()

# accession numbers seen in current but not prev
accessions = set()

# if True, only sequences present in the previous list will be output, instead of only seqs not present
output_existing_only = False

if len(sys.argv) >= 6:
    output_existing_only = (sys.argv[5].lower() == "true")

# read prev metadata
with open(sys.argv[2]) as f:
    for line_s in f:
        # just a list of accessions, not a JSON file
        prev_accessions.add(line_s.strip())

# read current metadata
with open(sys.argv[1]) as f:
    for line_s in f:
        line = json.loads(line_s)
        accession = get_if_exists(line, "accession")
        if accession == "" or (accession in prev_accessions and not output_existing_only) or (
                accession not in prev_accessions and output_existing_only):
            # skipping anything with no accession number (should not ever happen)
            # or that duplicates an old sequence (or that is new if output_existing_only is true)
            continue
        accessions.add(accession)
        loc = get_if_exists(line, "location", "geographicLocation")
        loc_div = ""
        if ":" in loc:
            loc_split = loc.split(":")
            loc = loc_split[0].strip()
            loc_div = loc_split[1].strip()
        output_lines[accession] = [
            accession,
            get_if_exists(line, "isolate", "name"),
            get_if_exists(line, "isolate", "collectionDate"),
            get_if_exists(line, "location", "geographicRegion"),
            loc,
            loc_div,
            get_if_exists(line, "releaseDate")
        ]

# read seqs
with open(sys.argv[3]) as f:
    with open(sys.argv[4], "w") as g:
        accession = None
        for line_s in f:
            line = line_s.split("\t")
            accession = line[0][1:].split()[0]
            if accession in accessions:
                g.write("\t".join(output_lines[accession] + [line[1].strip()]) + "\n")
