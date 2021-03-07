# frozen_string_literal: true

class DetailedGeoLocation < ApplicationRecord
  has_many :subscribed_geo_locations

  def self.unique_fields
    %i[region subregion division subdivision]
  end

  def cache_key(detailed_geo_location)
    DetailedGeoLocation.unique_fields.map { |f| detailed_geo_location.send(f) }.join
  end

  def self.existing_geo_location_ids_by_unique_fields
    @existing_geo_location_ids_by_unique_fields ||= Hash[DetailedGeoLocation.pluck([:id] + unique_fields)
                                                                            .map { |dg| [dg[1..], dg[0]] }.join]
  end

  def to_s
    "#{region}/#{subregion}/#{division}/#{subdivision}"
  end
end
