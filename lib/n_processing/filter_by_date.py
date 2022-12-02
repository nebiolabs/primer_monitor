import sys
import datetime
import logging
import re

date_fails = 0

if len(sys.argv) < 3:
    sys.stderr.write("usage: python filter_by_date.py <input TSV> <cutoff date>\n")
    sys.exit(1)

date_fields = [int(f) for f in sys.argv[2].split("-")]

cutoff_date = datetime.date(*date_fields)


def is_date_new(date_str):
    if date_str == "":
        return False
    split_date = date_str.split("-")
    rec_date = datetime.date(int(split_date[0]), int(split_date[1]), int(split_date[2]))
    return cutoff_date < rec_date


with open(sys.argv[1]) as f:
    for line_s in f:
        if line_s.startswith("ref_start"):
            continue
        line_s = line_s.strip()
        line = line_s.split("\t")
        if len(line) < 5:
            logging.error("unparseable line, skipping: " + line_s)
            continue
        valid_date = re.compile("\d{4}-\d{1,2}-\d{1,2}")
        if not valid_date.match(line[4].strip()):
            logging.error(
                "unparseable date, skipping: " + line[4] + "\nIn line: " + line_s
            )
            continue
        if line[4] == "":
            date_fails += 1
        if is_date_new(line[4]):
            print(line_s)
logging.warning("undated records: " + str(date_fails))
