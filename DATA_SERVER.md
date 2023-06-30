# igv.js Data Server Configuration

The igv.js visualization on the lineages page requires a data server from which to read the organism (currently SARS-CoV-2) 
reference sequence and BED files for primer sets and lineages. The URL of this data server is specified in the `.env`
file at the project root.

This server only needs to host static files and be able to accept HTTP Range requests, so something basic like NGINX 
should work for this purpose.

## Data Server File Paths

For the purposes of this section, we will leave off the domain of the IGV data server for clarity 
(e.g. a path below of `<organism name>/misc/tracks.json` indicates a full URL of `<IGV data server URL>/<organism name>/misc/tracks.json`).

Primer Monitor expects the data files hosted by the data server to be at the following paths:

### Metadata

* `<organism name>/misc/tracks.json`: A JSON file with key-value pairs mapping the internal name of a primer set (used in URLs),
to the name a user should see on the site (e.g. `"ARTIC_v4": "ARTIC v4"`). The list of available primer sets is
also obtained from the list of key-value pairs in this file.

* `<organism name>/misc/lineage_sets.json`: A JSON file with key-value pairs mapping the internal name of a lineage set 
(used in URLs), to the name a user should see on the site (e.g. `"BA.2": "BA.2.*"`). The list of available 
lineage sets is also obtained from the list of key-value pairs in this file.

### Primer sets

* `<organism name>/primer_sets/<primer set name>/<lineage group name>.bed`: A BED file (minimum 9 columns) containing a record for
each primer in the primer set. The color is used to mark whether the primer is affected by a variant in the specified
lineage group.

* `<organism name>/primer_sets_raw/<primer set name>.bed`: A BED file containing the raw primer annotations for
  each primer in the primer set.

### Lineage sets

* `<organism name>/lineage_sets/<lineage set>.txt`: A plain text file listing the Pango lineages in the database comprising
  a lineage set. The filenames match the internal names from `lineage_sets.json`.

* `<organism name>/lineage_variants/<lineage set>.bed`: A BED file listing the variants detected in the specified
  lineage set. The filenames match the internal names from `lineage_sets.json`.

* `<organism name>/lineage_variants/<lineage set>.bed.tbi`: A tabix index for `<lineage set>.bed`.

### Reference sequence

* `<organism name>/ref/<genbank accession>.fasta`: The reference sequence of the SARS-CoV-2 virus in FASTA format. For SARS-CoV-2, `<genbank accession>` is `NC_045512.2`.
* `<organism name>/ref/<genbank accession>.fasta.fai`: An index of the FASTA reference sequence.
* `<organism name>/ref/<genbank accession>.gff3`: A GFF3 file containing genome annotations.
* `<organism name>/ref/<genbank accession>.gff3.tbi`: A tabix index of the annotations file.

