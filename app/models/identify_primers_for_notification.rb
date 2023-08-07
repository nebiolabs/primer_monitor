# frozen_string_literal: true

class IdentifyPrimersForNotification < ApplicationRecord
  belongs_to :detailed_geo_location
  belongs_to :primer_set
  belongs_to :oligo

  self.primary_keys = :primer_set_id, :user_id, :oligo_id, :coords, :detailed_geo_location_id

  def coordinate
    coords
  end
end
