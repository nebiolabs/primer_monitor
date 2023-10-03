# frozen_string_literal: true

class GeoLocation < ApplicationRecord
  has_many :subscribed_geo_locations, dependent: :destroy

  def to_s
    "#{region}/#{division}"
  end
end
