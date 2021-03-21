# frozen_string_literal: true

class LocationAliasJoin < ApplicationRecord
  belongs_to :detailed_geo_location_alias
  belongs_to :detailed_geo_location
  has_many :fasta_records, foreign_key: :detailed_geo_location_id
end
