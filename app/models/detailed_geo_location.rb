# frozen_string_literal: true

class DetailedGeoLocation < ApplicationRecord
  has_many :fasta_records
  belongs_to :detailed_geo_location_alias

  def self.unique_fields
    %i[region subregion division subdivision]
  end

  def cache_key
    DetailedGeoLocation.unique_fields.map { |f| send(f) }.join('/')
  end

  def self.existing_geo_location_ids_by_unique_fields
    @existing_geo_location_ids_by_unique_fields ||=
      DetailedGeoLocation.pluck(:id, unique_fields.join(',')).each_with_object({}) do |dg_fields, h|
        h[dg_fields[1..].join('/')] = dg_fields[0]
      end
  end

  def to_s
    "#{region}/#{subregion}/#{division}/#{subdivision}"
  end
end
