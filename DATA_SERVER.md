# igv.js Data Server Configuration

The igv.js visualization on the lineages page requires a data server from which to read the SARS-CoV-2 
reference sequence and BED files for primer sets and lineages. The URL of this data server is specified in the `.env`
file at the project root.

This server only needs to host static files and be able to accept HTTP Range requests, so something basic like NGINX 
should work for this purpose.

## Data Server File Paths

For the purposes of this section, we will leave off the domain of the IGV data server for clarity 
(i.e. a path below of `/misc/tracks.json` indicates a full URL of `<IGV data server URL>/misc/tracks.json`).

Primer-monitor expects the data files hosted by the data server to be at the following paths:

### Metadata

* `/misc/tracks.json`: A JSON file with key-value pairs mapping the internal name of a primer set (used in URLs),
to the name a user should see on the site (e.g. `"ARTIC_v4": "ARTIC v4"`). The list of available primer sets is
also obtained from the list of key-value pairs in this file.

* `/misc/lineage_sets.json`: A JSON file with key-value pairs mapping the internal name of a lineage set 
(used in URLs), to the name a user should see on the site (e.g. `"BA.2": "BA.2.*"`). The list of available 
lineage sets is also obtained from the list of key-value pairs in this file.

### Primer sets

* `/primer_sets/<primer set name>/<lineage group name>.bed`: A BED file (minimum 9 columns) containing a record for
each primer in the primer set. The color is used to mark whether the primer is affected by a variant in the specified
lineage group.

### Reference sequence

* `/ref/NC_045512.2.fasta`: The reference sequence of the SARS-CoV-2 virus in FASTA format.
* `/ref/NC_045512.2.fasta.fai`: An index of the FASTA reference sequence.
* `/ref/NC_045512.2.gff3`: A GFF3 file containing SARS-CoV-2 genome annotations.

