# Finds all descendant lineages of the given designation
# Usage: curl <pangolin alias_key.json URL> | python get_lineage_names.py <lineage>

# this is also imported as a module by get_lineages_to_show.py

import json
import sys
import re
import datetime
import psycopg2
import os

def get_relevant_aliases(search_string, aliases):
    relevant_aliases = [search_string]
    for alias in aliases.keys():
        search_split = search_string.split(".")
        alias_split = alias.split(".")
        search_len = len(search_split)
        alias_len = len(alias_split)
        min_len = min(search_len, alias_len)
        if alias_split[:min_len] == search_split[:min_len]:
            extra_segment = ""
            # if search string is more specific than alias
            if search_len > alias_len:
                # add the extra specificity onto the alias
                extra_segment = ".".join([""]+search_split[alias_len:search_len])
            for child in aliases[alias]:
                relevant_aliases += get_relevant_aliases(child+extra_segment, aliases)
    return relevant_aliases

def process_aliases(alias_str, lineage_filename, search_string, conn):

    output = ""

    aliases_data = json.loads(alias_str)

    search_alias = re.sub('\.?\*$', '', search_string) # remove trailing ".*" or "*" from name

    reversed_aliases = {}

    for alias in aliases_data:
        if type(aliases_data[alias]) == type([]): # is a list
            for parent in aliases_data[alias]:
                if parent not in reversed_aliases:
                    reversed_aliases[parent] = []
                reversed_aliases[parent].append(alias)
        else:
            parent = aliases_data[alias]
            if parent not in reversed_aliases:
                reversed_aliases[parent] = []
            reversed_aliases[parent].append(alias)


    relevant_aliases = get_relevant_aliases(search_alias, reversed_aliases)

    relevant_aliases_for_search = []

    for relevant in relevant_aliases:
        relevant_aliases_for_search.append(relevant+".") # adding trailing "." so e.g. "A" doesn't match "AY"

    lineage_data = []
    with open(lineage_filename) as lineages_file:
        has_counts = False
        for line_s in lineages_file:
            if "," in line_s:
                # if the file has a second column
                line = line_s.split(",")
                if len(line) >= 2:
                    has_counts = True
            else:
                # make this a 1-element list so the rest of the code can stay the same
                line = [line_s.strip()]
            search_str = line[0]+"." # adding trailing "." for same reason as above
            for alias in relevant_aliases_for_search:
                #print(alias, search_str)
                if search_str.startswith(alias):
                    output+=(line[0]+"\n")
                    #print("checking "+str(line))
                    lineage_data.append(line[0])
                    break # break out of the inner loop and go to the next line
    num_seen = 0
    min_date = None
    max_date = None
    median_date = None
    if len(lineage_data) > 0 and conn is not None:
        cur = conn.cursor()
        cur.execute('SELECT COUNT(fasta_records.id) AS num_seen, min(fasta_records.date_submitted) AS first_seen, \
        max(fasta_records.date_submitted) AS last_seen, to_timestamp(percentile_cont(0.5) WITHIN GROUP \
        (ORDER BY cast(extract(epoch FROM date_submitted) AS integer)))::date AS median_date FROM fasta_records INNER JOIN pangolin_calls \
        ON pangolin_calls.id=fasta_records.pangolin_call_id INNER JOIN lineages ON lineages.id=pangolin_calls.lineage_id \
        WHERE lineages.name IN %s;', (tuple(lineage_data),))
        num_seen, min_date, max_date, median_date = cur.fetchone()
        #print(num_seen, min_date, max_date, median_date)
        cur.close()
    #print(str([output, num_seen, min_date, max_date, median_date]))
    return [output, num_seen, min_date, max_date, median_date]


if __name__ == "__main__":
    conn = psycopg2.connect("dbname="+os.environ['DB_NAME']+" user="+os.environ['DB_USER_RO']+" host="+os.environ['DB_HOST'])
    result = process_aliases(sys.stdin.read(), sys.argv[2], sys.argv[1], conn)
    conn.close()
    print(result[0])
    print("total,"+str(result[1]))
