if (($# < 7)); then
    echo "usage: ./qsub_pangolin.sh <working dir> <input file name> <pangolin path> <threads>" >&2;
    exit 1;
fi

$pwd=`pwd`;

qsub -S /bin/bash -pe smp $3 -N run_pangolin $pwd/run_pangolin.sh $pwd/$1 pangolin_calls.csv $pwd/pangolin_calls.done $pwd $2 $3 > $pwd/run_pangolin.log 2> $pwd/run_pangolin.err # qsub pangolin task

