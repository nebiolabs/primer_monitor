import sys
import math

NEAR_5P_SCORE = 1
NEAR_3P_SCORE = 2
MID_3P_SCORE = 4
OTHER_SCORE = 3
CUTOFF = int(sys.argv[4])
INDEL_MULTIPLIER = 10
LENGTH_EXPONENT_BASE = 2

primers_scored = {}

with open(sys.argv[1]) as f:
    for line_s in f:
        line = line_s.strip().split("\t")

        variant_start_ind = len(line) - 1
        primer_start_pos = int(line[1])
        primer_end_pos = int(line[2])
        strand = line[5]
        seq_name = line[3]

        variant_start_pos = None
        variant_end_pos = None
        variant_name = None

        while line[variant_start_ind] != line[0]:
            variant_start_ind -= 1
        try:
            variant_start_pos = int(line[variant_start_ind + 1])
            variant_end_pos = int(line[variant_start_ind + 2])
            variant_name = line[variant_start_ind + 3]
        except:
            print(line)
            raise

        cigar_op = variant_name.split("/")[0]
        overlap_start = (primer_start_pos if primer_start_pos > variant_start_pos else variant_start_pos)
        overlap_end = (primer_end_pos if primer_end_pos < variant_end_pos else variant_end_pos)
        overlap_len = overlap_end - overlap_start
        score = 0
        pos = overlap_start
        dist_5p = None
        dist_3p = None
        dir = (1 if strand == "+" else -1)
        if dir > 0:
            dist_5p = overlap_start - primer_start_pos
            dist_3p = primer_end_pos - overlap_start
        else:
            dist_5p = primer_end_pos - overlap_start
            dist_3p = overlap_start - primer_start_pos

        if cigar_op == "I" and overlap_start < primer_start_pos:  # insertion outside primer
            continue

        while pos <= overlap_end:
            new_score = 0
            if dist_5p <= 3:
                new_score = NEAR_5P_SCORE
            elif dist_3p <= 2:
                new_score = NEAR_3P_SCORE
            elif 2 < dist_3p <= 5:
                new_score = MID_3P_SCORE
            else:
                new_score = OTHER_SCORE

            dist_5p += dir
            dist_3p -= dir
            pos += 1
            if cigar_op == "D" or cigar_op == "I":  # indel
                new_score *= INDEL_MULTIPLIER

            score += new_score

        if seq_name not in primers_scored:
            primers_scored[seq_name] = line[:6]  # get BED record for primer
            primers_scored[seq_name][4] = 0  # init score

        primers_scored[seq_name][4] += (LENGTH_EXPONENT_BASE ** score)  # update score

    affected = open(sys.argv[2], "w")
    unaffected = open(sys.argv[3], "w")
    for primer in primers_scored:
        score = int(primers_scored[primer][4])
        ln_score = math.log(score)  # ln to reduce the range of scores
        if ln_score > 1000:
            ln_score = 1000  # cap at 1000 to follow BED spec
        primers_scored[primer][4] = str(ln_score)
        if score > CUTOFF:
            affected.write("\t".join(primers_scored[primer]) + "\n")
        else:
            unaffected.write("\t".join(primers_scored[primer]) + "\n")
    affected.close()
    unaffected.close()
