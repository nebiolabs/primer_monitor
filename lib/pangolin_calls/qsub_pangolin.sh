if (($# < 7)); then
    echo "usage: ./qsub_pangolin.sh <batch number> <working dir> <input file name> <output file name> <done file name> <pangolin path> <threads>" >&2;
    exit 1;
fi

qsub -S /bin/bash -pe smp $7 -N pango_batch$1 $2/run_pangolin.sh $2/$3 $4 $2/$5 $2 $6 $7 > $2/pango_batch$1.log 2> $2/pango_batch$1.err # qsub pangolin task


