# frozen_string_literal: true

class IdentifyPrimersForNotification < ApplicationRecord
  belongs_to :detailed_geo_location
  belongs_to :primer_set
  belongs_to :oligo

  def coordinate
    coords
  end
end
