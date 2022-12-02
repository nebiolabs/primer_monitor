import sys

if len(sys.argv) < 2:
    sys.stderr.write("usage: python cov_to_bedgraph.py <bedtools coverage file>\n")
    sys.exit(1)

with open(sys.argv[1]) as f:
    for line_s in f:
        line = line_s.split("\t")
        newline = "\t".join(
            [
                line[0],
                str(int(line[1]) - 1 + int(line[3])),
                str(int(line[1]) + int(line[3])),
                line[4],
            ]
        )
        print(newline.strip())
