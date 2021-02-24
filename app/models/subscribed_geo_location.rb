# frozen_string_literal: true

class SubscribedGeoLocation < ApplicationRecord
  belongs_to :user
  belongs_to :geo_location
  belongs_to :detailed_geo_location

  validates_uniqueness_of :user_id, scope: :geo_location_id
end
