# frozen_string_literal: true

class DetailedGeoLocation < ApplicationRecord
  has_many :subscribed_geo_locations

  def to_s
    "#{region}/#{subregion}/#{division}/#{subdivision}"
  end
end
