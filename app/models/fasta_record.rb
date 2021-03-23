# frozen_string_literal: true

class FastaRecord < ApplicationRecord
  belongs_to :detailed_geo_location

  def self.parse(metadata_tsv)
    raise "Unable to find counts file #{metadata_tsv}" unless File.exist?(metadata_tsv)

    new_fasta_records = []
    record_count = 0
    @new_locations = {}

    File.readlines(metadata_tsv).each do |line|
      next if line.start_with?("strain\t")

      record_count += 1
      record = build_fasta_record(line)
      new_fasta_records << record if record
    end
    @existing_fasta_strain_ids = nil # invalidates cache since any new records would not be present
    raise "Unable to parse any records from #{metadata_tsv}" if record_count.zero?

    new_fasta_records
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

    ActiveRecord::Base.logger.info("New fasta record: #{strain}")

    dg = get_dg(country, division, location, region)

    FastaRecord.new(strain: strain, gisaid_epi_isl: gisaid_epi_isl,
                    genbank_accession: genbank_accession, detailed_geo_location: dg,
                    date_collected: date, variant_name: variant_name)
  end

  def self.get_dg(country, division, location, region)
    # .presence converts each string to nil if it's empty
    new_dg = DetailedGeoLocation.new(world: 'World', region: region.presence, subregion: country.presence,
                                     division: division.presence, subdivision: location.presence)

    # fetches dg_ids from the cache if it already exists in the database, no harm if it's nil
    dg_id = DetailedGeoLocation.existing_geo_location_ids_by_unique_fields[new_dg.cache_key]
    dg = @new_locations[new_dg.cache_key] # re-use new locations

    if !dg_id && !dg
      # did not find an existing geolocation (dg_id) or a recently cached one (dg)
      new_dg.detailed_geo_location_alias = DetailedGeoLocationAlias.new_from_detailed_geolocation(new_dg)
      @new_locations[new_dg.cache_key] = new_dg # update new location with one that has an id
      ActiveRecord::Base.logger.info("New location: #{new_dg.cache_key}, id: #{new_dg.id}")
      dg = new_dg
    end

    dg
  end
end
