# frozen_string_literal: true

# Sequence counts placed on the map: each location's records count towards the map position of its state or country
# (see VariantPositions::CountryCentroids.point). Records in a location without a position are left out of the map.
module GeoPoints
  # {position name => {count:, point:}} for the records, counted by the location they were collected in
  def self.sums(records)
    counts = records.where.not(detailed_geo_location_id: nil).group(:detailed_geo_location_id).count
    locations = DetailedGeoLocation.where(id: counts.keys).index_by(&:id)
    counts.each_with_object({}) do |(location_id, count), sums|
      location = locations.fetch(location_id)
      point = VariantPositions::CountryCentroids.point(location.subregion, location.division)
      next if point.nil?

      sums[point[:name]] ||= { count: 0, point: }
      sums[point[:name]][:count] += count
    end
  end

  # one row per position with its count, the total sequences in it and the fraction, as the map draws them
  def self.rows(counts, totals)
    rows = counts.map do |name, entry|
      point = entry[:point]
      total = totals.dig(name, :count) || 0
      { name:, variant_count: entry[:count], sequence_count: total,
        fraction: total.zero? ? 0.0 : entry[:count].to_f / total, latitude: point[:lat], longitude: point[:lng] }
    end
    rows.sort_by { |row| [-row[:variant_count], row[:name]] }
  end
end
