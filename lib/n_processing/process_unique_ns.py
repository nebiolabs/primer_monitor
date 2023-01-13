import sys

region_n_counts = {}

seqs_to_dates = {}

regions_to_pos = {}

if len(sys.argv) < 2:
    sys.stderr.write(
        "usage: python process_unique_ns.py <sorted input BED file> [n fraction cutoff]\n"
    )
    sys.exit(1)

n_cutoff = 0

if len(sys.argv) >= 3:
    n_cutoff = float(sys.argv[2])

with open(sys.argv[1]) as f:
    for line in f:
        variant = [field.strip() for field in line.split("\t")]
        if variant[4] not in seqs_to_dates:
            seqs_to_dates[variant[4]] = variant[5]
        if variant[3] not in regions_to_pos:
            regions_to_pos[variant[3]] = [variant[1], variant[2]]
        if variant[3] not in region_n_counts:
            region_n_counts[variant[3]] = {variant[4]: int(variant[6])}
        elif variant[4] not in region_n_counts[variant[3]]:
            region_n_counts[variant[3]][variant[4]] = int(variant[6])
        else:
            region_n_counts[variant[3]][variant[4]] += int(variant[6])

for region in region_n_counts:
    for seq in region_n_counts[region]:
        n_frac = float(region_n_counts[region][seq]) / float(
            int(regions_to_pos[region][1]) - int(regions_to_pos[region][0])
        )
        if n_frac >= n_cutoff:
            region_num = region.split("_")[0]
            region_set_name = "_".join(region.split("_")[1:])
            print(
                "NC_045512.2\t"
                + regions_to_pos[region][0]
                + "\t"
                + regions_to_pos[region][1]
                + "\t"
                + region_num
                + "\t"
                + region_set_name
                + "\t"
                + seq
                + "\t"
                + seqs_to_dates[seq]
                + "\t"
                + str(region_n_counts[region][seq])
                + "\t"
                + str(n_frac)
            )
