# frozen_string_literal: true

module VariantPositions
  # Map position of a country, from config/geo/country_centroids.json. Matches the country names in the sequence
  # metadata, which vary in spelling, by their accent-free, lower-cased form.
  module CountryCentroids
    def self.normalize(name)
      I18n.transliterate(name.to_s).downcase.gsub(/[^a-z0-9]+/, ' ').strip
    end

    DATA = JSON.parse(Rails.root.join('config/geo/country_centroids.json').read)

    BY_NAME = DATA.fetch('countries').each_with_object({}) do |country, index|
      country.fetch('names').each { |name| index[normalize(name)] ||= country }
    end.freeze

    # the centroid entry ({'name', 'lat', 'lng', 'names'}) for a country name, or nil when it is not known
    def self.lookup(name)
      BY_NAME[normalize(name)]
    end
  end
end
