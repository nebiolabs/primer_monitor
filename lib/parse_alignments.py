
import sys
import re

def parse_cigar(cigar_string):
    parsed = re.findall("(\d+[\w=])", cigar_string)
    return parsed

# print("Sample\tReference position\tMismatch type\tMismatch")
for line in sys.stdin:
    if line.startswith("@"):
        continue
    line = line.strip().split("\t")
    query_name = line[0]
    ref_start = int(line[3])
    cigar = line[5]
    query_seq = line[9]
    parsed_cigar = parse_cigar(cigar)

    # print(parsed_cigar)
    query_start = 0
    read_consuming_ops = ("M", "I", "S", "=", "X")
    ref_consuming_ops = ("M", "D", "N", "=", "X")
    for item in parsed_cigar:
        length = int(item[:-1])
        type = item[len(item)-1:]
        if type == "D":
            output = "-" * length
        else:
            output = query_seq[query_start : query_start + length]
        if not type in "MHS=": #dont print match or clips
            print(f"{query_name}\t{ref_start}\t{type}\t{output}")
        if type in read_consuming_ops:
            query_start += length
        ref_start += length
        