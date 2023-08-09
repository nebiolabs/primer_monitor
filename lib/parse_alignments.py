"""
Parses variants from CIGAR strings in alignment SAM data.
"""
import sys
import re


def parse_cigar(cigar_string):
    """
    Parses CIGAR string into list of operations.
    :param cigar_string: CIGAR string
    :return: List of operations.
    """
    parsed = re.findall(r"(\d+[\w=])", cigar_string)
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
        op_type = item[len(item) - 1:]
        if op_type == "D":
            output = "-" * length
        else:
            output = query_seq[query_start: query_start + length]
        if op_type not in "MHS=":  # don't print match or clips
            print(f"{query_name}\t{ref_start}\t{op_type}\t{output}")
        if op_type in read_consuming_ops:
            query_start += length
        if op_type in ref_consuming_ops:
            ref_start += length
