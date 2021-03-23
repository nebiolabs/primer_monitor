# frozen_string_literal: true

class DetailedGeoLocationAlias < ApplicationRecord
  has_many :subscribed_geo_locations, inverse_of: :detailed_geo_location_alias
  has_many :detailed_geo_locations
  has_many :fasta_records, through: :detailed_geo_locations

  MIN_SEQUENCES_FOR_GEOLOCATION = 20
  scope :subscribable, lambda {
    joins(detailed_geo_locations: :fasta_records)
      .group('detailed_geo_location_aliases.id')
      .having("count(fasta_records.id) >= #{MIN_SEQUENCES_FOR_GEOLOCATION}")
      .order('world, region NULLS FIRST, subregion NULLS FIRST,
              subdivision NULLS FIRST, locality NULLS FIRST, sublocality NULLS FIRST')
  }

  def self.new_from_detailed_geolocation(detailed_geolocation)
    DetailedGeoLocationAlias.new(world: 'World', region: detailed_geolocation.region,
                                 subregion: detailed_geolocation.subregion,
                                 division: detailed_geolocation.division,
                                 subdivision: detailed_geolocation.subdivision)
  end

  def to_s
    "#{region}/#{subregion}/#{division}/#{subdivision}"
  end

  def name
    to_s
  end
end
