# frozen_string_literal: true

class FastaRecord < ApplicationRecord
  belongs_to :detailed_geo_location
  belongs_to :pangolin_call, foreign_key: :pangolin_call_id, optional: true

  def self.parse(metadata_tsv)
    raise "Unable to find counts file #{metadata_tsv}" unless File.exist?(metadata_tsv)

    new_fasta_records = []
    record_count = 0
    @new_locations = {}

    File.readlines(metadata_tsv).each do |line|
      next if line.start_with?("accession\t")

      record_count += 1
      record = build_fasta_record(line)
      new_fasta_records << record if record
    end
    @existing_fasta_accession_ids = nil # invalidates cache since any new records would not be present
    raise "Unable to parse any records from #{metadata_tsv}" if record_count.zero?

    new_fasta_records
  end

  def self.existing_fasta_accession_ids
    @existing_fasta_accession_ids ||= Hash[FastaRecord.pluck(:genbank_accession, :id)]
  end

  def self.build_fasta_record(line)
    (accession, strain, date, region, country, division, date_submitted) = line.chomp.split("\t")

    return unless strain && accession && region && country && division && date

    division = nil if division.blank?
    return if existing_fasta_accession_ids.key?(accession)

    ActiveRecord::Base.logger.info("New fasta record: \"#{accession}\"")
    # ensure that higher level locations also exist so users can select these for subscriptions
    find_or_create_dg_id(region, nil, nil, nil)
    find_or_create_dg_id(region, country, nil, nil)
    dg_id = find_or_create_dg_id(region, country, division, nil)

    FastaRecord.new(strain: strain,
                    genbank_accession: accession, detailed_geo_location_id: dg_id,
                    date_collected: date, date_submitted: date_submitted)
  end

  # fetches detailed geolocation record for the specified parameters, creates if necessary
  def self.find_or_create_dg_id(region, country, division, location)
    # .presence converts each string to nil if it's empty
    new_dg = DetailedGeoLocation.new(world: 'World', region: region.presence, subregion: country.presence,
                                     division: division.presence, subdivision: location.presence)

    # fetches dg_ids from the cache if it already exists in the database, no harm if it's nil
    dg_id = DetailedGeoLocation.existing_geo_location_ids_by_unique_fields[new_dg.cache_key]

    unless dg_id
      dg = @new_locations[new_dg.cache_key] # re-use new locations
      dg_id = dg.id if dg
    end

    unless dg_id
      # did not find an existing geolocation (dg_id)
      new_dg.detailed_geo_location_alias = DetailedGeoLocationAlias.new_from_detailed_geolocation(new_dg)
      new_dg.save!
      @new_locations[new_dg.cache_key] = new_dg # update new location with one that has an id
      ActiveRecord::Base.logger.info("New location: #{new_dg.cache_key}, id: #{new_dg.id}")
      dg_id = new_dg.id
    end

    dg_id
  end

end
