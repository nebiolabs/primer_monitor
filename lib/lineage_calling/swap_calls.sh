#!/usr/bin/env bash

# ensure this is relative to the script's location
source "$(dirname "$0")/../../.env";

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" >&2 -v "taxon_id=$1" <<CMDS;
DELETE FROM lineage_calls WHERE id IN (SELECT lineage_call_id FROM fasta_records
  WHERE fasta_records.lineage_call_id IS NOT NULL AND fasta_records.pending_lineage_call_id IS NOT NULL AND
   fasta_records.organism_taxon_id =(SELECT id FROM organism_taxa WHERE organism_taxa.ncbi_taxon_id=:'taxon_id'));
UPDATE fasta_records SET lineage_call_id=pending_lineage_call_id, pending_lineage_call_id=NULL
  WHERE fasta_records.pending_lineage_call_id IS NOT NULL;
CMDS

touch "$(basename "$PWD").complete_swap";