# frozen_string_literal: true

require 'set'

class LineageCall < ApplicationRecord
  belongs_to :lineage
  has_one :fasta_record, dependent: :nullify

  def self.parse(calls_csv, caller_id)
    raise "Unable to find calls file #{calls_csv}" unless File.exist?(calls_csv)

    new_calls = []
    record_count = 0

    File.readlines(calls_csv).each do |line|
      next if line.start_with?('taxon,')

      record_count += 1
      record = build_lineage_call(line, caller_id)
      new_calls << record if record
    end
    raise "Unable to parse any records from #{calls_csv}" if record_count.zero?

    new_calls
  end

  def self.build_lineage_call(line, caller_id)
    (taxon, lineage, metadata) = line.chomp.split(',', 3)

    return unless taxon && lineage

    lineage_rec = Lineage.find_by(name: lineage)

    raise "Failed to find lineage #{lineage}" if lineage_rec.nil?

    LineageCall.new(taxon: taxon, lineage_id: lineage_rec.id, lineage_callers_id: caller_id, metadata: metadata)
  end

  def self.update_fasta_recs(pending)
    id_field = pending ? 'pending_lineage_call_id' : 'lineage_call_id'

    # likely faster this way than with ActiveRecord
    # the interpolation is safe because it can only be one of the 2 hardcoded values above
    ActiveRecord::Base.connection.execute "UPDATE fasta_records SET #{id_field}=lineage_calls.id
  FROM lineage_calls WHERE lineage_calls.taxon=fasta_records.genbank_accession AND fasta_records.#{id_field} IS NULL;"
  end
end
