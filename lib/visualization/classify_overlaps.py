import sys
import math

NEAR_5P_SCORE = 1
NEAR_3P_SCORE = 2
MID_3P_SCORE = 4
OTHER_SCORE = 3
CUTOFF = int(sys.argv[3])
INDEL_MULTIPLIER = 10
LENGTH_EXPONENT_BASE = 2

primers_scored = {}

with open(sys.argv[1]) as f:
    for line_s in f:
        line = line_s.strip().split("\t")
        variant_start_ind = len(line)-1
        variant_name = None
        while line[variant_start_ind] != line[0]:
            variant_start_ind-=1
        try:
            line[1] = int(line[1])
            line[2] = int(line[2])
            line[variant_start_ind+1] = int(line[variant_start_ind+1])
            line[variant_start_ind+2] = int(line[variant_start_ind+2])
            variant_name = line[variant_start_ind+3]
        except:
            print(line)
            raise
        cigar_op = variant_name.split("/")[0]
        overlap_start = (line[1] if line[1] > line[variant_start_ind+1] else line[variant_start_ind+1])
        overlap_end = (line[2] if line[variant_start_ind+2] < line[variant_start_ind+2] else line[variant_start_ind+2])
        overlap_len = overlap_end-overlap_start
        score = 0
        pos = overlap_start
        dist_5p = None
        dist_3p = None
        dir = (1 if line[5] == "+" else -1)
        new_score = 0
        if dir > 0:
            dist_5p = overlap_start-line[1]
            dist_3p = line[2]-overlap_start
        else:
            dist_5p = line[2]-overlap_start
            dist_3p = overlap_start-line[1]

        if cigar_op == "I" and overlap_start<int(line[1]): # insertion outside primer
            continue

        while pos <= overlap_end:
            if dist_5p <= 3:
                new_score = NEAR_5P_SCORE
            elif dist_3p <= 2:
                new_score = NEAR_3P_SCORE
            elif dist_3p >2 and dist_3p <= 5:
                new_score = MID_3P_SCORE
            else:
                new_score = OTHER_SCORE

            dist_5p+=dir
            dist_3p-=dir
            pos+=1
            if cigar_op == "D" or cigar_op == "I": # indel
                new_score = new_score * INDEL_MULTIPLIER

            score += new_score

        if line[3] not in primers_scored:
            primers_scored[line[3]]=line[:6]
            primers_scored[line[3]][4] = 0

        primers_scored[line[3]][4]+=(LENGTH_EXPONENT_BASE**score)

    affected = open(sys.argv[2], "w")
    for primer in primers_scored:
        score = primers_scored[primer][4]
        if score > CUTOFF:
            primers_scored[primer][4] = math.log(primers_scored[primer][4]) # ln to reduce the range of scores
            if primers_scored[primer][4] > 1000:
                score = 1000 # cap at 1000 to follow BED spec
            i = 0
            for field in primers_scored[primer]: # convert non-string fields to strings before .join("\t")
                if type(primers_scored[primer][i]) != type("a string"):
                    primers_scored[primer][i] = str(primers_scored[primer][i])
                i+=1
            affected.write("\t".join(primers_scored[primer])+"\n")
    affected.close()