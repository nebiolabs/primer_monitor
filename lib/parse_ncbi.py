import sys
import json

# returns the value from within the dict, or "" if a key is missing
def get_if_exists(cur_dict, *keys):
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

# read prev metadata
with open(sys.argv[2]) as f:
    for line_s in f:
        line = json.loads(line_s)
        accession = get_if_exists(line, "accession")
        if accession == "":
            continue # skipping anything with no accession number (should not ever happen)
        prev_accessions.add(accession)

# read current metadata
with open(sys.argv[1]) as f:
    for line_s in f:
        line = json.loads(line_s)
        accession = get_if_exists(line, "accession")
        if accession == "" or accession in prev_accessions:
            continue # skipping anything with no accession number (should not ever happen) or that duplicates an old sequence
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
            get_if_exists(line, "releaseDate"),
            "" # sequence
            ]

# read seqs
with open(sys.argv[3]) as f:
    with open(sys.argv[4], "w") as g:
        g.write("accession\tstrain\tdateCollected\tregion\tcountry\tdivision\tdateReleased\tsequence\n")

        accession = None
        for line_s in f:
            line = line_s.split("\t")
            accession = line[0][1:].split()[0]
            if accession in accessions:
                g.write("\t".join(output_lines[accession].append(line[1].strip())+"\n"))
