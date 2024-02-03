"""
Gets "interesting" lineages to display on the igv.js page
"interesting" = "has at least <cutoff> seqs in it"

Usage: get_lineages_to_show.py <root lineages> <lineages CSV> <seq counts CSV> <lineage sets path> [overrides path]

root lineages: a comma-separated string of root lineage names from which to start the tree traversal (e.g. A,B)

lineages CSV: the CSV of lineage data from get_lineage_data.sh

seq counts CSV: a file of day-by-day sequence counts for each lineage

lineage sets path: the path where lineage set definition files are output

overrides path: a path to a file listing a set of Pango lineage names to always consider "interesting"
"""
import argparse
import get_lineage_names
import sys
import datetime
import psycopg2
import os
import logging

MIN_SEQS = 50


def seqs_in_date_range(seqs_per_day, start, end):
    """

    Returns the total number of new sequences in a given date range.

    :param seqs_per_day: Data on how many new sequences were added each day.
    :param start: The start date of the range.
    :param end: The end date of the range.
    :return: The total sequence count.
    """
    total = 0
    num_days = (end - start).days
    for i in range(num_days):
        the_date = (start + datetime.timedelta(days=i)).strftime('%Y-%m-%d')
        if the_date in seqs_per_day:
            total += seqs_per_day[the_date]
    return total


def passes_filter(lineage, total, seqs_per_day, lineage_metadata):
    """

    Does this lineage qualify as "interesting"?

    :param lineage: The lineage group name.
    :param total: The total sequence count of this lineage group.
    :param seqs_per_day: Data on the number of new seqs each day.
    :param lineage_metadata: Lineage metadata
    :return: True if the lineage group is "interesting", otherwise False
    """

    today = datetime.datetime.now().date()
    if lineage_metadata[lineage]['last_seen'] is None:
        # no recorded lineages, fail
        return False

    # get day-by-day seqs submitted (date_submitted field) and use that data to compute based on fractions

    # total number of sequences (for any lineage) during the relevant timeframe
    total_seqs = seqs_in_date_range(seqs_per_day, lineage_metadata[lineage]['first_seen'],
                                    lineage_metadata[lineage]['last_seen'])

    seqs_after_median = seqs_in_date_range(seqs_per_day, lineage_metadata[lineage]['median_date'],
                                           lineage_metadata[lineage]['last_seen'])

    # increase cutoff for long-lasting lineages so they don't overwhelm new ones
    base_cutoff = 0.05 * total_seqs
    base_cutoff_after_median = 0.05 * seqs_after_median

    time_diff_median = (today - lineage_metadata[lineage]['median_date']).days
    time_diff_newest = (today - lineage_metadata[lineage]['last_seen']).days

    old_adjustment = 1
    # massive score penalty for really old lineages to get rid of stuff like B.1.*
    if (time_diff_newest > 90 and time_diff_median > 180) or time_diff_median > 360:
        old_adjustment = 100 + (max(time_diff_newest - 90, 0))

    logging.debug("testing full: " + lineage + ", seen: " + str(total) + ", overall: " + str(total_seqs)
                  + ", cutoff:" + str(base_cutoff) + ", from:" + str(lineage_data[lineage]['first_seen'])
                  + ", to:" + str(lineage_data[lineage]['last_seen']))

    logging.debug("testing after median: " + lineage + ", seen: " + str(total / 2)
                  + ", overall: " + str(seqs_after_median) + ", cutoff:" + str(base_cutoff_after_median)
                  + ", from:" + str(lineage_data[lineage]['median_date']) + ", to:" + str(
        lineage_data[lineage]['last_seen']))

    full_range = total > base_cutoff * old_adjustment and total >= MIN_SEQS
    after_median_range = total / 2 > base_cutoff_after_median * old_adjustment and total >= MIN_SEQS

    logging.debug("result: " + str(full_range or after_median_range))

    return full_range or after_median_range


def get_descendant_lineages(result, base):
    """
    Gets all descendant lineages for a given lineage listed in the lineages file,
    pulled from a defined set fo lineages.
    :param result: The set of lineages to search through.
    :param base: The base lineage to look for descendants of.
    :return: A list of lineage "prefixes" for descendants.
    """
    base_depth = len(base.split("."))
    prefixes = set()
    # loop over lineage names
    for rec in result.split("\n"):
        line = rec.strip()
        # if this is a child of the starting lineage
        if line.startswith(base + "."):
            # add 1 level deeper than base
            prefix = ".".join(line.split(".")[:base_depth + 1])
            prefixes.add(prefix)
        # if this is something else (e.g. an alias)
        elif not base.startswith(line):
            # add the prefix
            prefix = line.split(".")[0]
            prefixes.add(prefix)
        else:
            # don't add anything
            pass
    return prefixes


def get_lineages(lineage_tree, lineages_file, starting_lineages, seqs_per_day, lineage_metadata, conn, lineages=None,
                 root_lineages=None):
    """
    Recursively traverses the lineages tree to find "interesting" lineages.

    :param lineage_tree: The RSV lineage data as a dictionary
    :param lineages_file: A file listing all lineages in the database
    :param starting_lineages: The "root" lineages to start from for this recursive call
    :param seqs_per_day: Data on the number of new sequences each day.
    :param lineage_metadata: Lineage data pulled from the database.
    :param conn: Database connection.
    :param lineages: The resulting list of lineages, accumulated here with each recursive call.
    :param root_lineages: The list of *original* root lineages, to pass it to subsequent recursive calls.
    :return: A list of potentially-interesting lineages.
    """
    # if this is the initial call
    if lineages is None:
        lineages = set()
    if root_lineages is None:
        # set the root_lineages list so neither base clade is dropped
        root_lineages = set(starting_lineages)
    for lineage in starting_lineages:
        # get total seqs in this lineage group, and all lineage names
        result, total, min_date, max_date, median_date = get_lineage_names.process_aliases(lineage_tree, lineages_file,
                                                                                           lineage, conn, organism_slug)
        # if this lineage is "interesting" (has more than cutoff # of seqs or is a root lineage)
        if lineage not in lineage_metadata:
            lineage_metadata[lineage] = {'num_seen': total,
                                         'first_seen': min_date,
                                         'last_seen': max_date,
                                         'median_date': median_date}
        if (passes_filter(lineage, total, seqs_per_day,
                          lineage_metadata) or lineage in root_lineages) and lineage not in lineages:
            # add it and repeat recursively on that lineage
            lineages.add(lineage)
            lineages = get_lineages(lineage_tree, lineages_file, get_descendant_lineages(result, lineage), seqs_per_day,
                                    lineage_metadata, conn, lineages, root_lineages)
    return lineages


parser = argparse.ArgumentParser()
parser.add_argument('root_lineages')
parser.add_argument('organism_slug')
parser.add_argument('lineages_csv')
parser.add_argument('seq_counts_csv')
parser.add_argument('output_path')
parser.add_argument('overrides_path', nargs='?')
parsed_args = parser.parse_args()

base_lineages = parsed_args.root_lineages.split(",")
lineage_list_file = parsed_args.lineages_csv
seqs_per_day_file = parsed_args.seq_counts_csv
organism_slug = parsed_args.organism_slug

rsv_lineage_data = {}

cur_lineage = None
for line in sys.stdin:
    if line.startswith("## "):
        cur_lineage = line.split(" ")[-1].strip()
    if " * parent:" in line:
        parent = line.split(" ")[-1][1:].split("]")[0].strip()
        if parent == 'none':
            continue
        lineage_split = cur_lineage.split(".")
        parent_split = parent.split(".")
        # does the child lineage simply have extra elements
        min_len = min(len(parent_split), len(lineage_split))
        if lineage_split[:min_len] != parent_split[:min_len]:
            # get rid of common suffix to get an alias
            # e.g. this will convert "B.D.4.1.1 child of B.D.E.1"
            # to "B.D.4.1 alias of B.D.E"
            while lineage_split[-1] == parent_split[-1]:
                lineage_split.pop()
                parent_split.pop()
            alias_cur = ".".join(lineage_split)
            alias_parent = ".".join(parent_split)
            if cur_lineage in rsv_lineage_data:
                rsv_lineage_data[cur_lineage].append(parent)
            else:
                rsv_lineage_data[cur_lineage] = [parent]

daily_new_seq_counts = {}

with open(seqs_per_day_file) as f:
    for line_s in f:
        daily_seq_rec = line_s.strip().split(",")
        daily_new_seq_counts[daily_seq_rec[0]] = int(daily_seq_rec[1])

output_path = parsed_args.output_path

overrides = []

if parsed_args.overrides_path is not None:
    overrides_file = parsed_args.overrides_path
    with open(overrides_file) as f:
        for override_line in f:
            overrides.append(override_line.strip())

lineage_data = {}
with open(parsed_args.lineages_csv) as f:
    for line_s in f:
        lineage_rec = line_s.strip().split(",")
        lineage_data[lineage_rec[0]] = {'num_seen': lineage_rec[1],
                                        'first_seen': datetime.datetime.strptime(lineage_rec[2], '%Y-%m-%d').date(),
                                        'last_seen': datetime.datetime.strptime(lineage_rec[3], '%Y-%m-%d').date(),
                                        'median_date': datetime.datetime.strptime(lineage_rec[4], '%Y-%m-%d').date()}

db_conn = psycopg2.connect(
    "dbname=" + os.environ['DB_NAME'] + " user=" + os.environ['DB_USER_RO'] + " host=" + os.environ['DB_HOST'])
interesting_lineages = get_lineages(rsv_lineage_data, lineage_list_file, base_lineages, daily_new_seq_counts,
                                    lineage_data, db_conn)
db_conn.close()

for override_lineage in overrides:
    interesting_lineages.add(override_lineage)

lineage_groups = {}

reversed_aliases = get_lineage_names.reverse_aliases(rsv_lineage_data)

print("{")
for lineage_group in interesting_lineages:
    recorded_lineages = get_lineage_names.process_aliases(rsv_lineage_data, lineage_list_file, lineage_group, None,
                                                          organism_slug)[0]
    # if this lineage group has a single alias name, use that name
    display_name = lineage_group + "*"
    if lineage_group in reversed_aliases and len(reversed_aliases[lineage_group]) == 1:
        display_name = reversed_aliases[lineage_group][0] + "* (" + lineage_group + "*)"
    with open(output_path + "/" + lineage_group + ".txt", "w") as f:
        f.write(recorded_lineages)
    print('"' + lineage_group + '": "' + display_name + '",')
print('"all": "All"')
print("}")
