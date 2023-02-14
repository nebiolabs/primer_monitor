#!/usr/bin/awk -f
BEGIN {
    all = 0
    if(ARGV[1] == "all")
    {
        all = 1
        ARGV[1] = "" # skip this "file"
    }

    strainsfile = ARGV[2]
    ARGV[2] = "-" # read data from stdin
    ARGC = 3
}

FILENAME == ARGV[1] && all != 1 { #lineages data to read
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
