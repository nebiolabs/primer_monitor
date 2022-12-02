mkfifo $$_unique_regions;
mkfifo $$_unique_ns;
python3 get_unique_regions.py $1 $2 > $$_unique_regions & # get the regions that belong to only one amplicon
bedtools intersect -wo -a $$_unique_regions -b $3 | cut -f 1,2,3,4,8,9,12 > $$_unique_ns & # get Ns in unique regions
python process_unique_ns.py $$_unique_ns $4; # process this to get the number and fraction of Ns in the region and handle the N fraction cutoff
rm $$_unique_regions $$_unique_ns;