# frozen_string_literal: true

class SubscribedGeoLocation < ApplicationRecord
  belongs_to :user
  belongs_to :detailed_geo_location_alias

  validates_uniqueness_of :user_id, scope: :detailed_geo_location_alias_id

  def to_s
    detailed_geo_location_alias.name
  end
end
