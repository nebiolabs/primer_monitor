"""
Gets "interesting" lineages to display on the igv.js page
"interesting" = "has at least <cutoff> seqs in it"
"""
import get_lineage_names
import sys
import json
import datetime
import psycopg2
import os

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
    # time_range = (((lineage_data[lineage]['last_seen']-lineage_data[lineage]['first_seen']).days)/90)
    base_cutoff = 0.05 * total_seqs  # cutoff*max(time_range, 0.5)
    base_cutoff_after_median = 0.05 * seqs_after_median

    time_diff_median = (today - lineage_metadata[lineage]['median_date']).days
    time_diff_newest = (today - lineage_metadata[lineage]['last_seen']).days

    old_adjustment = 1
    # massive score penalty for really old lineages to get rid of stuff like B.1.*
    if (time_diff_newest > 90 and time_diff_median > 180) or time_diff_median > 360:
        old_adjustment = 100 + (max(time_diff_newest - 90, 0))

    # print("testing full: "+lineage+", seen: "+str(total)+", overall: "+str(total_seqs)\
    # +", cutoff:"+str(base_cutoff)+", from:"+str(lineage_data[lineage]['first_seen'])\
    # +", to:"+str(lineage_data[lineage]['last_seen']))

    # print("testing after median: "+lineage+", seen: "+str(total/2)\
    # +", overall: "+str(seqs_after_median)+", cutoff:"+str(base_cutoff_after_median)
    # +", from:"+str(lineage_data[lineage]['median_date'])+", to:"+str(lineage_data[lineage]['last_seen']))

    full_range = total > base_cutoff * old_adjustment and total >= MIN_SEQS
    after_median_range = total / 2 > base_cutoff_after_median * old_adjustment and total >= MIN_SEQS

    # print("result: "+str(full_range or after_median_range))

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


def get_lineages(alias_str, lineages_file, starting_lineages, seqs_per_day, lineage_metadata, conn, lineages=None,
                 root_lineages=None):
    """
    Recursively traverses the lineages tree to find "interesting" lineages.

    :param alias_str: The Pango aliases_key file as a string
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
        result, total, min_date, max_date, median_date = get_lineage_names.process_aliases(alias_str, lineages_file,
                                                                                           lineage, conn)
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
            lineages = get_lineages(alias_str, lineages_file, get_descendant_lineages(result, lineage), seqs_per_day,
                                    lineage_metadata, conn, lineages, root_lineages)
    return lineages


base_lineages = sys.argv[1].split(",")
pango_aliases_data = sys.stdin.read()
lineage_list_file = sys.argv[2]
seqs_per_day_file = sys.argv[3]

daily_new_seq_counts = {}

with open(seqs_per_day_file) as f:
    for line_s in f:
        daily_seq_rec = line_s.strip().split(",")
        daily_new_seq_counts[daily_seq_rec[0]] = int(daily_seq_rec[1])

output_path = sys.argv[4]
lineage_data = {}
with open(sys.argv[2]) as f:
    for line_s in f:
        lineage_rec = line_s.strip().split(",")
        lineage_data[lineage_rec[0]] = {'num_seen': lineage_rec[1],
                                        'first_seen': datetime.datetime.strptime(lineage_rec[2], '%Y-%m-%d').date(),
                                        'last_seen': datetime.datetime.strptime(lineage_rec[3], '%Y-%m-%d').date(),
                                        'median_date': datetime.datetime.strptime(lineage_rec[4], '%Y-%m-%d').date()}

db_conn = psycopg2.connect(
    "dbname=" + os.environ['DB_NAME'] + " user=" + os.environ['DB_USER_RO'] + " host=" + os.environ['DB_HOST'])
interesting_lineages = get_lineages(pango_aliases_data, lineage_list_file, base_lineages, daily_new_seq_counts,
                                    lineage_data, db_conn)
db_conn.close()

lineage_groups = {}

aliases_data = json.loads(pango_aliases_data)
reversed_aliases = get_lineage_names.reverse_aliases(aliases_data)

print("{")
for lineage_group in interesting_lineages:
    recorded_lineages = get_lineage_names.process_aliases(pango_aliases_data, lineage_list_file, lineage_group, None)[0]
    # if this lineage group has a single alias name, and that name is not a recombinant lineage (starts with X), use that name
    if lineage_group in reversed_aliases and len(reversed_aliases[lineage_group]) == 1 and not reversed_aliases[lineage_group][0].startswith("X"):
        lineage_group = reversed_aliases[lineage_group][0]
    with open(output_path + "/" + lineage_group + ".txt", "w") as f:
        f.write(recorded_lineages)
    print('"' + lineage_group + '*": "' + lineage_group + '",')
print("}")
