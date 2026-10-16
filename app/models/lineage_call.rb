# frozen_string_literal: true

class LineageCall < ApplicationRecord
  belongs_to :lineage
  has_one :fasta_record, dependent: :nullify

  def self.parse(calls_csv, caller_id)
    raise "Unable to find calls file #{calls_csv}" unless File.exist?(calls_csv)

    new_calls = []
    record_count = 0

    File.readlines(calls_csv).each do |line|
      next if line.start_with?('taxon,', 'seqName,')

      record_count += 1
      record = build_lineage_call(line, caller_id)
      new_calls << record if record
    end
    raise "Unable to parse any records from #{calls_csv}" if record_count.zero?

    new_calls
  end

  def self.build_lineage_call(line, caller_id)
    (taxon, lineage, metadata) = line.chomp.split(',', 4)

    return unless taxon && lineage

    lineage_rec = Lineage.find_by(name: lineage)

    raise "Failed to find lineage #{lineage}" if lineage_rec.nil?

    LineageCall.new(taxon:, lineage_id: lineage_rec.id, lineage_caller_id: caller_id, metadata:)
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
