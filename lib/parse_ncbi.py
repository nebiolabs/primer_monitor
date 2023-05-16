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

        do_print = False
        accession = None
        for line in f:
            if line.startswith(">"): # FASTA header
                if accession is not None and do_print:
                    g.write("\t".join(output_lines[accession])+"\n")
                if line[1:].split()[0] in accessions:
                    do_print = True
                    accession = line[1:].split()[0]
                else:
                    do_print = False
            else:
                if do_print:
                    output_lines[accession][-1] += line.strip()

        if accession is not None:
            g.write("\t".join(output_lines[accession])+"\n")