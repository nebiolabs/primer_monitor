if (($# < 2)); then
    echo "usage: ./qsub_pangolin.sh <input file name> <threads>" >&2;
    exit 1;
fi

workdir=`pwd`;
pangolin $workdir/$1 -t $2 -o $workdir --outfile pangolin_calls.csv; # run pangolin
touch $workdir/pangolin_calls.done; # create file to indicate this is done