BEGIN {
    all = 0
    if(ARGV[1] == "all")
    {
        all = 1
    }
    for(i=1; i<ARGC; i++)
    {
        lineages[ARGV[i]] = 1 #just defining the key so I can use "in"
        delete ARGV[i] # getting it out of ARGV so it doesn't try to read it as a file
    }
}

$5 in lineages || all == 1 {
   print
}