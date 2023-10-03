"""
Filters duplicate sequences out of an NCBI dataset given a list of sequences in the database.
Usage: ./parse_ncbi.py <metadata> <old accessions> <new sequences> <output filename>
"""
import json
import argparse


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


parser = argparse.ArgumentParser()
parser.add_argument("new_metadata", help="The current day's NCBI metadata file.")
parser.add_argument("new_sequences", help="The current day's NCBI sequences file.")
parser.add_argument("old_accessions", help="A file with the list of accessions already in the database.")
parser.add_argument("output_filename", help="The name of the output file.")

# if True, only sequences present in the previous list will be output, instead of only seqs not present
parser.add_argument("-e", "--output-existing",
                    help="Outputs sequences already in the database instead of new sequences", action="store_true")

args = parser.parse_args()

# maps accession numbers to data
output_lines = {}

# accession numbers seen in prev
prev_accessions = set()

# accession numbers seen in current but not prev
accessions = set()

# read prev metadata
# the prev metadata is just a list of accessions as it's output from a DB query
with open(args.old_accessions) as f:
    for line_s in f:
        prev_accessions.add(line_s.strip())

# read current metadata
with open(args.new_metadata) as f:
    for line_s in f:
        line = json.loads(line_s)
        accession = get_if_exists(line, "accession")
        if accession == "" or (accession in prev_accessions and not args.output_existing) or (
                accession not in prev_accessions and args.output_existing):
            # skipping anything with no accession number (should not ever happen)
            # or that duplicates an old sequence (or that is new if args.output_existing is true)
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
with open(args.new_sequences) as f:
    with open(args.output_filename, "w") as g:
        accession = None
        for line_s in f:
            line = line_s.split("\t")
            accession = line[0][1:].split()[0]
            if accession in accessions:
                g.write("\t".join(output_lines[accession] + [line[1].strip()]) + "\n")
