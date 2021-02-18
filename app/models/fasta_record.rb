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
      _paper_url, _date_submitted, variant_name) = line.chomp.split("\t")

    return unless strain && gisaid_epi_isl && genbank_accession && region && country && division && date && variant_name
    return if FastaRecord.exists?(strain: strain)
    region = strain.split("/")[0]
    division = strain.split("/")[1].split("-")[0]

    # The geo location record needs to exist before the fasta record does
    unless GeoLocation.exists?(region: region, division: division)
      GeoLocation.new(region: region, division: division).save!
    end

    geo_location_id = GeoLocation.find_by(region: region, division: division).id

    FastaRecord.new(strain: strain, gisaid_epi_isl: gisaid_epi_isl,
                    genbank_accession: genbank_accession, geo_location_id: geo_location_id,
                    date_collected: date, variant_name: variant_name)
  end
end
