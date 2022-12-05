import sys
import subprocess
import os

if len(sys.argv) < 4: # missing args
    sys.stderr.write("usage: python batch_seqs.py <batch size> <pangolin path> <pangolin threads>")
    sys.exit(1)

groupsize = int(sys.argv[1])/2 # divided by 2 because a seq is 2 lines

procs = [] # list of procs
running = [] # is this process running
filenames = [] # base filename for proc
outputs = [] # output data from proc

running_count = 0 # how many procs are running

i = 0 # how many seqs in current batch
line_group="" # current batch
groupnum = 0 # number of current batch

# a pangolin header, so they can be removed from outputs other than the first
header_line = "taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,scorpio_notes,version,pangolin_version,scorpio_version,constellation_version,is_designated,qc_status,qc_notes,note\n"

def run_group(): # runs the current batch
    global running_count
    # base file name for this batch
    filename = "tmp_"+str(os.getpid())+"_group_"+str(groupnum)
    with open(filename+".txt", "w") as g:
        g.write(line_group) # save seqs to file
    cwd = os.getcwd()
    # run pangolin
    p = subprocess.Popen(['bash', 'qsub_pangolin.sh', str(groupnum), cwd, filename+".txt", filename+"_out.txt", filename+"_done.txt", sys.argv[2], sys.argv[3]])
    # add process to lists
    procs.append(p)
    running.append(True)
    running_count+=1
    outputs.append(None)
    filenames.append(filename)

for line in sys.stdin:
    if i==groupsize: 
        # if there is a full batch of seqs, run it
        run_group()
        # reset counters
        i=0
        groupnum+=1
        # start new batch
        line_group = ""
    line_group+=line
    i+=1

 # if there is anything left over, run it
if i>0:
    run_group()

while running_count > 0: # repeatedly check until all procs finished
    for i in range(len(procs)):
         # if process is flagged as running but its done file exists (i.e. it just finished)
        if running[i] and os.path.exists(filenames[i]+"_done.txt"):
            running[i] = False
            running_count-=1
            with open(filenames[i]+"_out.txt") as f:
                outputs[i] = f.read()  # read output
                if i>0: # remove headers from all except the 1st file
                    outputs[i] = outputs[i][len(header_line):]
            # remove temp files
            os.remove(filenames[i]+".txt")
            os.remove(filenames[i]+"_out.txt")
            os.remove(filenames[i]+"_done.txt")

# combine outputs
for output in outputs:
    print(output, end="")
