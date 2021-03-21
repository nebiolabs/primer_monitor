# frozen_string_literal: true

class DetailedGeoLocationAlias < ApplicationRecord
  has_many :subscribed_geo_locations
  has_many :detailed_geo_locations
  has_many :fasta_records, through: :detailed_geo_locations

  def self.new_from_detailed_geolocation(detailed_geolocation)
    DetailedGeoLocationAlias.new(world: 'World', region: detailed_geolocation.region,
                                 subregion: detailed_geolocation.subregion,
                                 division: detailed_geolocation.division,
                                 subdivision: detailed_geolocation.subdivision)
  end
end
