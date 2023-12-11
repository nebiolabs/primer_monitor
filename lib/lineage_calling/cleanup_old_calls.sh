#!/usr/bin/env bash

# ensure this is relative to the script's location
source "$(dirname "$0")/../../.env";

# deletes any lineage call not connected to a fasta record

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" >&2 -v "taxon_id=$1" <<CMDS;
DELETE FROM lineage_calls WHERE NOT EXISTS
  (SELECT 1 FROM fasta_records WHERE
  fasta_records.organism_taxon_id =(SELECT id FROM organism_taxa WHERE organism_taxa.ncbi_taxon_id=:'taxon_id')
  AND fasta_records.lineage_call_id=lineage_calls.id
  OR fasta_records.pending_lineage_call_id=lineage_calls.id);
CMDS