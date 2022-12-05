import sys
import subprocess
import os

if len(sys.argv) < 4:
    sys.stderr.write("usage: python batch_seqs.py <batch size> <pangolin path> <pangolin threads>")
    sys.exit(1)

groupsize = int(sys.argv[1])/2 # divided by 2 because a seq is 2 lines

procs = []

running = []

filenames = []

outputs = []

running_count = 0

i = 0
line_group=""
groupnum = 0

header_line = "taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,scorpio_notes,version,pangolin_version,scorpio_version,constellation_version,is_designated,qc_status,qc_notes,note\n"

def run_group():
    global running_count
    filename = "tmp_"+str(os.getpid())+"_group_"+str(groupnum)
    with open(filename+".txt", "w") as g:
        g.write(line_group)
    cwd = os.getcwd()
    p = subprocess.Popen(['bash', 'qsub_pangolin.sh', str(groupnum), cwd, filename+".txt", filename+"_out.txt", filename+"_done.txt", sys.argv[2], sys.argv[3]])
    procs.append(p)
    running.append(True)
    running_count+=1
    outputs.append(None)
    filenames.append(filename)

for line in sys.stdin:
    if i==groupsize:
        run_group()
        i=0
        groupnum+=1
        line_group = ""
    line_group+=line
    i+=1

if i>0:
    run_group()

while running_count > 0:
    for i in range(len(procs)):
        if running[i] and os.path.exists(filenames[i]+"_done.txt"):
            running[i] = False
            running_count-=1
            with open(filenames[i]+"_out.txt") as f:
                outputs[i] = f.read()
                if i>0:
                    outputs[i] = outputs[i][len(header_line):]
            os.remove(filenames[i]+".txt")
            os.remove(filenames[i]+"_out.txt")
            os.remove(filenames[i]+"_done.txt")

for output in outputs:
    print(output, end="")
