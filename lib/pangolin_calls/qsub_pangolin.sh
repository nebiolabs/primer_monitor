qsub -S /bin/bash -pe smp $7 -N pango_batch$1 $2/run_pangolin.sh $2/$3 $4 $2/$5 $2 $6 $7 > $2/pango_batch$1.log 2> $2/pango_batch$1.err


