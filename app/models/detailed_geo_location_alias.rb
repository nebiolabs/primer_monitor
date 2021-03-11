# frozen_string_literal: true

class DetailedGeoLocation < ApplicationRecord
    has_many :subscribed_geo_locations
    has_one :location_alias_joins
end