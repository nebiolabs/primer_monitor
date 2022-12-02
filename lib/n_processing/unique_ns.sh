mkfifo $$_variants_new; # make named pipe
if (($# == 5)); then # if N variants file provided
    cat $5 > $$_variants_new & # load it in
elif (($# == 4)); then # otherwise, query for N variants
    ./extract_ns.sh $4 > $$_variants_new &
else # if anything else, error and print usage
    echo "usage: unique_ns.sh <primer BED> <primer set name> <N fraction cutoff> <date cutoff (ignored if file given)> [N variants file]" >&2;
    rm $$_variants_new; # clean up named pipe before exiting
    exit 1;
fi

# make more named pipes
mkfifo $$_unique_regions;
mkfifo $$_unique_ns;

python3 get_unique_regions.py $1 $2 > $$_unique_regions & # get the regions that belong to only one amplicon
bedtools intersect -wo -a $$_unique_regions -b $$_variants_new | cut -f 1,2,3,4,8,9,12 > $$_unique_ns & # get Ns in unique regions
python process_unique_ns.py $$_unique_ns $3; # process this to get the number and fraction of Ns in the region and handle the N fraction cutoff
rm $$_variants_new $$_unique_regions $$_unique_ns; # clean up the named pipes
