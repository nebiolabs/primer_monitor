# frozen_string_literal: true

class DetailedGeoLocationAlias < ApplicationRecord
  has_many :subscribed_geo_locations, inverse_of: :detailed_geo_location_alias, dependent: :destroy
  has_many :detailed_geo_locations, dependent: :destroy
  has_many :fasta_records, through: :detailed_geo_locations

  MIN_SEQUENCES_FOR_GEOLOCATION = 20
  scope :with_enough_sequences, lambda {
    joins(detailed_geo_locations: :fasta_records)
      .group('detailed_geo_location_aliases.id')
      .having("count(fasta_records.id) >= #{MIN_SEQUENCES_FOR_GEOLOCATION}")
      .order('world, region NULLS FIRST, subregion NULLS FIRST, division NULLS FIRST,
              subdivision NULLS FIRST, locality NULLS FIRST, sublocality NULLS FIRST')
  }
  # we are interested in locations that are named, but where there there is no
  # precision below these as well  e.g. United States (not but a specific state)
  scope :world, -> { where(region: nil) }
  scope :regions, -> { where(subregion: nil).where.not(region: nil) }
  scope :subregions, -> { where(division: nil).where.not(subregion: nil) }

  def self.subscribable
    DetailedGeoLocationAlias.world +
      DetailedGeoLocationAlias.regions +
      DetailedGeoLocationAlias.subregions +
      DetailedGeoLocationAlias.with_enough_sequences
  end

  def self.new_from_detailed_geolocation(detailed_geolocation)
    DetailedGeoLocationAlias.new(world: 'World', region: detailed_geolocation.region,
                                 subregion: detailed_geolocation.subregion,
                                 division: detailed_geolocation.division,
                                 subdivision: detailed_geolocation.subdivision)
  end

  def to_s
    if region.nil?
      world
    else
      [region, subregion, division, subdivision].compact.join('/')
    end
  end

  def name
    to_s
  end
end
