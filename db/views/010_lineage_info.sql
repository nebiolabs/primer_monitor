
CREATE MATERIALIZED VIEW public.lineage_info AS
    SELECT lineages.name,
          lineages.organism_id,
          count(fasta_records.id) AS times_seen,
          COALESCE(max(fasta_records.date_collected), max(fasta_records.date_submitted)) AS last_seen,
          COALESCE(min(fasta_records.date_collected), min(fasta_records.date_submitted)) AS first_seen
   FROM lineages
            JOIN lineage_calls ON lineage_calls.lineage_id = lineages.id
            JOIN fasta_records ON fasta_records.lineage_call_id = lineage_calls.id
   GROUP BY lineages.id
WITH DATA;