# Gets "interesting" lineages to display on the igvjs page
# "interesting" = "has at least <cutoff> seqs in it"

import get_lineage_names
import sys
import json
import datetime
import psycopg2
import os

MIN_SEQS = 50

def seqs_in_date_range(seqs_per_day, start, end):
    total = 0
    num_days = (end-start).days
    for i in range(num_days):
        the_date = (start+datetime.timedelta(days=i)).strftime('%Y-%m-%d')
        if the_date in seqs_per_day:
            total += seqs_per_day[the_date]
    return total


def passes_filter(lineage, total, seqs_per_day, lineage_data):
    today = datetime.datetime.now().date()
    if lineage_data[lineage]['last_seen'] is None:
        # no recorded lineages, fail
        return False

    # get day-by-day seqs submitted (date_submitted field) and use that data to compute based on fractions

    # total number of sequences (for any lineage) during the relevant timeframe
    total_seqs = seqs_in_date_range(seqs_per_day, lineage_data[lineage]['first_seen'], lineage_data[lineage]['last_seen'])

    seqs_after_median = seqs_in_date_range(seqs_per_day, lineage_data[lineage]['median_date'], lineage_data[lineage]['last_seen'])

    # increase cutoff for long-lasting lineages so they don't overwhelm new ones
    #time_range = (((lineage_data[lineage]['last_seen']-lineage_data[lineage]['first_seen']).days)/90)
    base_cutoff = 0.05*total_seqs #cutoff*max(time_range, 0.5)
    base_cutoff_after_median = 0.05*seqs_after_median


    time_diff_median = (today-lineage_data[lineage]['median_date']).days
    time_diff_newest = (today-lineage_data[lineage]['last_seen']).days

    old_adjustment = 1
    # massive score penalty for really old lineages to get rid of stuff like B.1.*
    if (time_diff_newest > 90 and time_diff_median > 180) or time_diff_median > 360:
        old_adjustment = 100+(max(time_diff_newest-90, 0))

    #print("testing full: "+lineage+", seen: "+str(total)+", overall: "+str(total_seqs)+", cutoff:"+str(base_cutoff)+", from:"+str(lineage_data[lineage]['first_seen'])+", to:"+str(lineage_data[lineage]['last_seen']))

    #print("testing after median: "+lineage+", seen: "+str(total/2)+", overall: "+str(seqs_after_median)+", cutoff:"+str(base_cutoff_after_median)+", from:"+str(lineage_data[lineage]['median_date'])+", to:"+str(lineage_data[lineage]['last_seen']))

    full_range = total > base_cutoff*old_adjustment and total >= MIN_SEQS
    after_median_range = total/2 > base_cutoff_after_median*old_adjustment and total >= MIN_SEQS

    #print("result: "+str(full_range or after_median_range))

    return full_range or after_median_range

# gets all sublineages for a given lineage listed in the lineages file
def get_sublineages(result, base):
    base_depth = len(base.split("."))
    prefixes = set()
    # loop over lineage names
    for rec in result.split("\n"):
        line = rec.strip()
        # if this is a child of the starting lineage
        if line.startswith(base+"."):
            # add 1 level deeper than base
            prefix = ".".join(line.split(".")[:base_depth+1])
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

# recursively traverses the lineages tree to find "interesting" lineages
def get_lineages(alias_str, lineages_file, starting_lineages, seqs_per_day, lineage_data, conn, lineages = set(), root_lineages = None):
    # if this is the initial call
    if root_lineages is None:
        # set the root_lineages list so neither base clade is dropped
        root_lineages = set(starting_lineages)
    for lineage in starting_lineages:
        # get total seqs in this lineage group, and all lineage names
        result, total, min_date, max_date, median_date = get_lineage_names.process_aliases(alias_str, lineages_file, lineage, conn)
        # if this lineage is "interesting" (has more than cutoff # of seqs or is a root lineage)
        if lineage not in lineage_data:
            lineage_data[lineage] = {'num_seen':total,
                                     'first_seen':min_date,
                                     'last_seen':max_date,
                                     'median_date':median_date}
        if (passes_filter(lineage, total, seqs_per_day, lineage_data) or lineage in root_lineages) and lineage not in lineages:
            # add it and repeat recursively on that lineage
            lineages.add(lineage)
            lineages = get_lineages(alias_str, lineages_file, get_sublineages(result, lineage), seqs_per_day, lineage_data, conn, lineages, root_lineages)
    return lineages

starting_lineages = sys.argv[1].split(",")
alias_str = sys.stdin.read()
lineages_file = sys.argv[2]
seqs_per_day_file = sys.argv[3]

seqs_per_day = {}

with open(seqs_per_day_file) as f:
    for line_s in f:
        line = line_s.strip().split(",")
        seqs_per_day[line[0]] = int(line[1])

output_path = sys.argv[4]
lineage_data = {}
with open(sys.argv[2]) as f:
    for line_s in f:
        line = line_s.strip().split(",")
        lineage_data[line[0]] = {'num_seen':line[1],
                                 'first_seen':datetime.datetime.strptime(line[2], '%Y-%m-%d').date(),
                                 'last_seen':datetime.datetime.strptime(line[3], '%Y-%m-%d').date(),
                                 'median_date':datetime.datetime.strptime(line[4], '%Y-%m-%d').date()}

conn = psycopg2.connect("dbname="+os.environ['DB_NAME']+" user="+os.environ['DB_USER_RO']+" host="+os.environ['DB_HOST'])
interesting_lineages = get_lineages(alias_str, lineages_file, starting_lineages, seqs_per_day, lineage_data, conn)
conn.close()

lineage_groups = {}

aliases_data = json.loads(alias_str)
reversed_aliases = reverse_aliases(aliases_data)

print("{")
for lineage in interesting_lineages:
    lineages = get_lineage_names.process_aliases(alias_str, lineages_file, lineage, None)[0]
    if lineage in reversed_aliases:
        lineage = reversed_aliases[lineage]
    with open(output_path+"/"+lineage+".txt", "w") as f:
        f.write(lineages)
    print('"'+lineage+'*": "'+lineage+'",')
print("}")

