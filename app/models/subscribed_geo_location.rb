# frozen_string_literal: true

class SubscribedGeoLocation < ApplicationRecord
  belongs_to :user
  belongs_to :detailed_geo_location_alias
  has_many :proposed_notifications, dependent: :destroy

  validates :user_id, uniqueness: { scope: :detailed_geo_location_alias_id }

  def to_s
    detailed_geo_location_alias.name
  end
end
