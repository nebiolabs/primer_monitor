# frozen_string_literal: true

class LocationAliasJoin < ApplicationRecord
    belongs_to :DetailedGeoLocationAlias
    belongs_to :DetailedGeoLocation
end
