BEGIN {
    FS="\t"
    OFS="\t"
    filename_metadata = cur_date ".metadata"
    filename_seqs = cur_date ".fasta"
}

$1 != "accession" {
    print ">" $1 "\n" $8 > filename_seqs;
}

{
    print $1, $2, $3, $4, $5, $6, $7 > filename_metadata;
}