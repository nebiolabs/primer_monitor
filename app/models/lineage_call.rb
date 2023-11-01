# frozen_string_literal: true

require 'set'

class LineageCall < ApplicationRecord
  belongs_to :lineage
  has_one :fasta_record, dependent: :nullify

  def self.parse(pangolin_csv)
    raise "Unable to find calls file #{pangolin_csv}" unless File.exist?(pangolin_csv)

    new_calls = []
    record_count = 0

    File.readlines(pangolin_csv).each do |line|
      next if line.start_with?('taxon,')

      record_count += 1
      record = build_pangolin_call(line)
      new_calls << record if record
    end
    raise "Unable to parse any records from #{pangolin_csv}" if record_count.zero?

    new_calls
  end

  def self.build_pangolin_call(line)
    (taxon, lineage, conflict, ambiguity_score, scorpio_call, scorpio_support, scorpio_conflict,
      scorpio_notes, version, pangolin_version, scorpio_version, constellation_version,
      is_designated, qc_status, qc_notes, note) = line.chomp.split(',')

    unless taxon && lineage && version && pangolin_version && scorpio_version && constellation_version && is_designated
      return
    end

    lineage_rec = Lineage.find_by(name: lineage)

    raise "Failed to find lineage #{lineage}" if lineage_rec.nil?

    LineageCall.new(taxon: taxon, lineage_id: lineage_rec.id, conflict: conflict,
                    ambiguity_score: ambiguity_score, scorpio_call: scorpio_call,
                    scorpio_support: scorpio_support, scorpio_conflict: scorpio_conflict,
                    scorpio_notes: scorpio_notes, version: version, pangolin_version: pangolin_version,
                    scorpio_version: scorpio_version, constellation_version: constellation_version,
                    is_designated: is_designated,
                    qc_status: qc_status, qc_notes: qc_notes, note: note)
  end

  def self.update_fasta_recs(pending)
    id_field = pending ? 'pending_pangolin_call_id' : 'pangolin_call_id'

    # likely faster this way than with ActiveRecord
    # the interpolation is safe because it can only be one of the 2 hardcoded values above
    ActiveRecord::Base.connection.execute "UPDATE fasta_records SET #{id_field}=pangolin_calls.id
  FROM pangolin_calls WHERE pangolin_calls.taxon=fasta_records.genbank_accession AND fasta_records.#{id_field} IS NULL;"
  end
end
