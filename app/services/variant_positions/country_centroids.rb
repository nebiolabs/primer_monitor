# frozen_string_literal: true

module VariantPositions
  # Map positions for sequence locations, from config/geo: a country's centroid, or for the countries with state data
  # (currently the United States) the state's centroid. Country names in the sequence metadata vary in spelling, so they
  # are matched by their accent-free, lower-cased form.
  module CountryCentroids
    def self.normalize(name)
      I18n.transliterate(name.to_s).downcase.gsub(/[^a-z0-9]+/, ' ').strip
    end

    COUNTRIES = JSON.parse(Rails.root.join('config/geo/country_centroids.json').read).fetch('countries')
    STATES = JSON.parse(Rails.root.join('config/geo/us_state_centroids.json').read)

    BY_NAME = COUNTRIES.each_with_object({}) do |country, index|
      country.fetch('names').each { |name| index[normalize(name)] ||= country }
    end.freeze

    STATE_COUNTRY = normalize(STATES.fetch('country'))
    STATES_BY_NAME = STATES.fetch('states').index_by { |state| normalize(state.fetch('name')) }.freeze

    # the country's centroid entry ({'name', 'lat', 'lng', 'names'}) for a country name, or nil when it is not known
    def self.lookup(name)
      BY_NAME[normalize(name)]
    end

    # where a location appears on the map: {name:, lat:, lng:} for its state when the country has state data and the
    # division is a state, otherwise for its country; nil when the country is not known
    def self.point(subregion, division)
      country = lookup(subregion)
      return if country.nil?

      state = STATES_BY_NAME[normalize(division)] if normalize(country['name']) == STATE_COUNTRY
      if state
        { name: "#{state['name']}, #{country['name']}", lat: state['lat'], lng: state['lng'] }
      else
        { name: country['name'], lat: country['lat'], lng: country['lng'] }
      end
    end
  end
end
