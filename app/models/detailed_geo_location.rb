# frozen_string_literal: true

class GeoLocation < ApplicationRecord
    has_many :subscribed_geo_locations
  
    def to_s
      "#{region}/#{subregion}/#{division}/#{subdivision}/#{locality}/#{latitude}/#{longitude}"
    end
  end
  