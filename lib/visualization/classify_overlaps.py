import sys
import math

NEAR_5P_SCORE = 1
NEAR_3P_SCORE = 2
MID_3P_SCORE = 4
OTHER_SCORE = 3
CUTOFF = int(sys.argv[3])
DELETION_MULTIPLIER = 10
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
        fivep_dist = None
        threep_dist = None
        dir = (1 if line[5] == "+" else -1)
        new_score = 0
        if dir > 0:
            fivep_dist = overlap_start-line[1]
            threep_dist = line[2]-overlap_start
        else:
            fivep_dist = line[2]-overlap_start
            threep_dist = overlap_start-line[1]

        if cigar_op == "I" and overlap_start<int(line[1]): # insertion outside primer
            continue

        while pos <= overlap_end:
            if fivep_dist <= 3:
                new_score = NEAR_5P_SCORE
            elif threep_dist <= 2:
                new_score = NEAR_3P_SCORE
            elif threep_dist >2 and threep_dist <= 5:
                new_score = MID_3P_SCORE
            else:
                new_score = OTHER_SCORE

            fivep_dist+=dir
            threep_dist-=dir
            pos+=1
            if cigar_op == "D": # deletion
                new_score = new_score * DELETION_MULTIPLIER

            score += new_score

        if line[3] not in primers_scored:
            primers_scored[line[3]]=line[:6]
            primers_scored[line[3]][4] = 0

        primers_scored[line[3]][4]+=(LENGTH_EXPONENT_BASE**score)

    affected = open(sys.argv[2], "w")
    for primer in primers_scored:
        score = primers_scored[primer][4]
        if score > CUTOFF:
            print(str(score)+" > "+str(CUTOFF))
            i = 0
            for field in primers_scored[primer]:
                primers_scored[primer][i] = str(primers_scored[primer][i])
                i+=1
            primers_scored[primer][4] = math.log(primers_scored[primer][4]) # log to reduce the range of scores
            if primers_scored[primer][4] > 1000:
                score = 1000 # cap at 1000 to follow BED spec
            print("primers_scored[primer][4]  = " + primers_scored[primer][4])
            affected.write("\t".join(primers_scored[primer])+"\n")
    affected.close()