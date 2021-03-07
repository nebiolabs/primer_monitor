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
    @existing_fasta_strain_ids = nil # invalidates cache since any new records would not be present
    raise "Unable to parse any records from #{metadata_tsv}" if record_count.zero?

    metadata
  end

  def self.existing_fasta_strain_ids
    @existing_fasta_strain_ids ||= Hash[FastaRecord.pluck(:strain, :id)]
  end

  def self.build_fasta_record(line)
    (strain, _virus, gisaid_epi_isl, genbank_accession, date, region, country, division, location,
      _region_exposure, _country_exposure, _division_exposure, _segment, _length, _host, _age, _sex,
      _pangolin_lineage, _gisaid_clade, _originating_lab, _submitting_lab, _authors, _url, _title,
      _paper_url, _date_submitted, variant_name) = line.chomp.split("\t")

    return unless strain && gisaid_epi_isl && genbank_accession && region && country && division && location && date
    return if existing_fasta_strain_ids.key?(strain)

    # Converts each string to nil if it's empty
    region = region.presence
    country = country.presence
    division = division.presence
    location = location.presence

    dg = DetailedGeoLocation.new(world: 'World', region: region, subregion: country, division: division,
                                 subdivision: location)

    # fetches dg_id from the cache if it already exists in the database, no harm if it's nil
    dg_id = DetailedGeoLocation.existing_geo_location_ids_by_unique_fields[DetailedGeoLocation.cache_key(dg)]

    fa = FastaRecord.new(strain: strain, gisaid_epi_isl: gisaid_epi_isl,
                         genbank_accession: genbank_accession, detailed_geo_location_id: dg_id,
                         date_collected: date, variant_name: variant_name)

    # dg_id was not found in the database, add the new record so it will be saved during import
    fa.detailed_geo_location = dg unless dg_id
  end
end
