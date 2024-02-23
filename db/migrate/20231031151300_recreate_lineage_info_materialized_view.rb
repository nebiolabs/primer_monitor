class RecreateLineageInfoMaterializedView < ActiveRecord::Migration[6.1]
  def up
    execute <<-SQL
      CREATE MATERIALIZED VIEW lineage_info AS SELECT lineages.name, lineages.organism_id, COUNT(fasta_records.id) AS times_seen,
      COALESCE(MAX(fasta_records.date_collected), MAX(fasta_records.date_submitted)) as last_seen,
      COALESCE(MIN(fasta_records.date_collected), MIN(fasta_records.date_submitted)) as first_seen
      FROM lineages INNER JOIN lineage_calls ON lineage_calls.lineage_id=lineages.id
      INNER JOIN fasta_records ON fasta_records.lineage_call_id=lineage_calls.id GROUP BY lineages.id;
      GRANT SELECT ON lineage_info TO primer_monitor_ro;
    SQL
  end

  def down
    execute <<-SQL
      DROP MATERIALIZED VIEW lineage_info;
    SQL
  end
end
