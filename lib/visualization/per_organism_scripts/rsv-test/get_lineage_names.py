"""
Finds all descendant lineages of the given designation
Usage: curl <pangolin alias_key.json URL> | python get_lineage_names.py <lineage>

this is also imported as a module by get_lineages_to_show.py
"""

import json
import sys
import re
import psycopg2
import os
import logging


def get_relevant_aliases(search_string, aliases):
    """
    Gets all "unaliased" names for a lineage.
    """
    relevant_aliases = [search_string]
    search_split = search_string.split(".")
    for alias in aliases.keys():
        alias_split = alias.split(".")
        search_len = len(search_split)
        alias_len = len(alias_split)
        min_len = min(search_len, alias_len)
        if alias_split[:min_len] == search_split[:min_len]:
            extra_segment = ""
            # if search string is more specific than alias
            if search_len > alias_len:
                # add the extra specificity onto the alias
                extra_segment = ".".join([""] + search_split[alias_len:search_len])
            for child in aliases[alias]:
                relevant_aliases += get_relevant_aliases(child + extra_segment, aliases)
    return relevant_aliases


def reverse_aliases(aliases_data):
    """
    Reverses an aliases JSON from alias->original to be original->alias
    """
    reversed_aliases = {}

    for alias in aliases_data:
        if isinstance(aliases_data[alias], list):  # is a list
            for parent in aliases_data[alias]:
                if parent not in reversed_aliases:
                    reversed_aliases[parent] = []
                reversed_aliases[parent].append(alias)
        else:
            parent = aliases_data[alias]
            if parent not in reversed_aliases:
                reversed_aliases[parent] = []
            reversed_aliases[parent].append(alias)

    return reversed_aliases


def process_aliases(lineage_parents, lineage_filename, search_string, db_conn):
    """
    Recursively finds all descendant lineages of a given lineage recorded in the database.
    """
    output = ""

    reversed_aliases = reverse_aliases(lineage_parents)

    search_alias = re.sub(r'\.?\*$', '', search_string)  # remove trailing ".*" or "*" from name

    relevant_aliases = get_relevant_aliases(search_alias, reversed_aliases)

    relevant_aliases_for_search = []

    for relevant in relevant_aliases:
        relevant_aliases_for_search.append(relevant + ".")  # adding trailing "." so e.g. "A" doesn't match "AY"

    lineage_data = []
    with open(lineage_filename) as lineages_file:
        for line_s in lineages_file:
            if "," in line_s:
                # if the file has a second column
                line = line_s.split(",")
            else:
                # make this a 1-element list so the rest of the code can stay the same
                line = [line_s.strip()]
            search_str = line[0] + "."  # adding trailing "." for same reason as above
            for alias in relevant_aliases_for_search:
                logging.debug(alias, search_str)
                if search_str.startswith(alias):
                    output += (line[0] + "\n")
                    logging.debug("checking "+str(line))
                    lineage_data.append(line[0])
                    break  # break out of the inner loop and go to the next line
    num_seen = 0
    min_date = None
    max_date = None
    median_date = None
    if len(lineage_data) > 0 and db_conn is not None:
        cur = db_conn.cursor()
        cur.execute('SELECT COUNT(fasta_records.id) AS num_seen, min(fasta_records.date_submitted) AS first_seen, \
        max(fasta_records.date_submitted) AS last_seen, to_timestamp(percentile_cont(0.5) WITHIN GROUP \
        (ORDER BY cast(extract(epoch FROM date_submitted) AS integer)))::date AS median_date \
        FROM fasta_records INNER JOIN lineage_calls ON lineage_calls.id=fasta_records.lineage_call_id \
        INNER JOIN lineages ON lineages.id=lineage_calls.lineage_id \
        WHERE lineages.name IN %s;', (tuple(lineage_data),))
        num_seen, min_date, max_date, median_date = cur.fetchone()
        logging.debug(num_seen, min_date, max_date, median_date)
        cur.close()
    logging.debug(str([output, num_seen, min_date, max_date, median_date]))
    return [output, num_seen, min_date, max_date, median_date]


if __name__ == "__main__":
    conn = psycopg2.connect(
        "dbname=" + os.environ['DB_NAME'] + " user=" + os.environ['DB_USER_RO'] + " host=" + os.environ['DB_HOST'])
    result = process_aliases(sys.stdin.read(), sys.argv[2], sys.argv[1], conn)
    conn.close()
    print(result[0])
    print("total," + str(result[1]))
