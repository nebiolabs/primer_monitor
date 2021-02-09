# frozen_string_literal: true

class FastaRecord < ApplicationRecord
  def self.parse(metadata_tsv)
    raise "Unable to find counts file #{metadata_tsv}" unless File.exist?(metadata_tsv)

    metadata = []
    record_count = 0

    File.readlines(metadata_tsv).each do |line|
      record_count += 1
      record = build_fasta_record(line)
      metadata << record if record
    end

    raise "Unable to parse any records from #{metadata_tsv}" if record_count.zero?

    metadata
  end

  def self.build_fasta_record(line)
    (strain, _virus, gisaid_epi_isl, genbank_accession, date, region, country, division, _location,
      _region_exposure, _country_exposure, _division_exposure, _segment, _length, _host, _age, _sex,
      _pangolin_lineage, _gisaid_clade, _originating_lab, _submitting_lab, _authors, _url, _title,
      _paper_url, _date_submitted) = line.chomp.split("\t")

    return unless strain && gisaid_epi_isl && genbank_accession && region && country && division && date
    return if FastaRecord.exists?(strain: strain)
    region = strain.split("/")[0]
    division = strain.split("/")[1]

    FastaRecord.new(strain: strain, gisaid_epi_isl: gisaid_epi_isl,
                    genbank_accession: genbank_accession, region: region,
                    country: country, division: division,
                    date_submitted: date)
  end
end
