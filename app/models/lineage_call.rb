# frozen_string_literal: true

class LineageCall < ApplicationRecord
  belongs_to :lineage
  has_one :fasta_record, dependent: :nullify

  # Returns unsaved LineageCall records for each row of calls_csv. Every lineage named in the file must
  # already exist; Lineage.parse and its import run first.
  def self.parse(calls_csv, caller_id)
    # Resolve every lineage name in one query rather than one round trip per row.
    lineage_ids = Lineage.ids_by_name(Lineage.lineage_names(calls_csv))

    new_calls = []
    File.foreach(calls_csv) do |line|
      next if line.start_with?('taxon,', 'seqName,')

      record = build_lineage_call(line, caller_id, lineage_ids)
      new_calls << record if record
    end

    new_calls
  end

  def self.build_lineage_call(line, caller_id, lineage_ids)
    (taxon, lineage, metadata) = line.chomp.split(',', 4)

    return unless taxon && lineage

    lineage_id = lineage_ids[lineage]

    raise "Failed to find lineage #{lineage}" if lineage_id.nil?

    LineageCall.new(taxon:, lineage_id:, lineage_caller_id: caller_id, metadata:)
  end

  # Points fasta records at the calls just imported (call_ids). Limiting this to the batch keeps the parallel
  # load_lineage_data tasks of a recall off each other's rows. A recall (pending) re-calls every record, so it
  # overwrites any pending call a failed earlier recall left behind; nightly calls only fill records without one.
  def self.update_fasta_recs(pending, call_ids)
    return if call_ids.empty?

    id_field = pending ? 'pending_lineage_call_id' : 'lineage_call_id'
    only_unset = pending ? '' : "AND fasta_records.#{id_field} IS NULL"

    # likely faster this way than with ActiveRecord
    # the interpolation is safe because it can only be one of the hardcoded values above
    ActiveRecord::Base.connection.execute sanitize_sql_array([<<~SQL.squish, call_ids])
      UPDATE fasta_records SET #{id_field}=lineage_calls.id
      FROM lineage_calls WHERE lineage_calls.id IN (?)
      AND lineage_calls.taxon=fasta_records.genbank_accession #{only_unset}
    SQL
  end
end
