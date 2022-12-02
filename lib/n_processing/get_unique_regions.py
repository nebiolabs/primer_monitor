import sys
import re

regions = []

unique_regions = []

primers = []

last_primer_start = -1

amplicons = {}

amplicon_nums = []

if len(sys.argv) < 3:
    sys.stderr.write(
        "usage: python get_unique_regions.py <sorted input BED file> <primer set name>\n"
    )
    sys.exit(1)

set_name = sys.argv[2]

with open(sys.argv[1]) as f:
    regex = re.compile("(\d+)_(LEFT|RIGHT)")
    for line in f:
        primer = [field.strip() for field in line.split("\t")]
        match = regex.search(primer[3])
        if match is None:
            sys.stderr.write("fatal: malformed primer name '" + primer[3] + "'\n")
            sys.exit(1)
        match_groups = match.groups()
        if match_groups[0] not in amplicon_nums:
            amplicon_nums.append(match_groups[0])
        if match_groups[0] not in amplicons:
            amplicons[match_groups[0]] = [[], []]
        if match_groups[1] == "LEFT":
            amplicons[match_groups[0]][0].append(int(primer[1]))
        else:
            amplicons[match_groups[0]][1].append(int(primer[2]))

i = 0
for amplicon_id in amplicon_nums:
    amplicon = amplicons[amplicon_id]
    if len(amplicon_nums) == 1:
        new_amplicon = [max(amplicon[0]), min(amplicon[1])]
    elif i == 0:
        new_amplicon = [
            max(amplicon[0]),
            min(min(amplicons[amplicon_nums[i + 1]][0]), min(amplicon[1])),
            amplicon_id + "_" + set_name,
        ]
    elif i == len(amplicon_nums) - 1:
        new_amplicon = [
            max(max(amplicons[amplicon_nums[i - 1]][1]), max(amplicon[0])),
            min(amplicon[1]),
            amplicon_id + "_" + set_name,
        ]
    else:
        new_amplicon = [
            max(max(amplicons[amplicon_nums[i - 1]][1]), max(amplicon[0])),
            min(min(amplicons[amplicon_nums[i + 1]][0]), min(amplicon[1])),
            amplicon_id + "_" + set_name,
        ]
    i += 1
    unique_regions.append(new_amplicon)

for region in unique_regions:
    print("NC_045512.2\t" + ("\t".join([str(field) for field in region])))
