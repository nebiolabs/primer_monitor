#!/usr/bin/awk -f

# filter sequences by lineage

# Usage: cat <variants bed file> | filter_variants.awk <lineage set file> <sequences file>
# Variants bed file must be passed on stdin (shown above with cat for clarity)

BEGIN {
    all = 0
    if(ARGV[1] == "all")
    {
        all = 1
        ARGV[1] = "" # skip this "file"
    }

    strainsfile = ARGV[2] # this argument is a file to write to, not read from
    ARGV[2] = "-" # read data from stdin
    ARGC = 3 # adjust argument count to account for removed argument
}

FILENAME == ARGV[1] && all != 1 { #lineage_variants data to read
    lineages[$0] = 1
}

FILENAME == ARGV[2] { # BED records, FILENAME is "-"
    if ($5 in lineages || all == 1) # either in lineage list or set to "all"
    {
        print
        strains[$4] = 1 # record strain seen
    }
}

END {
    for (strain in strains)
    {
        print strain >> strainsfile
    }
}
