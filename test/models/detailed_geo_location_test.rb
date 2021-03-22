# frozen_string_literal: true

require 'test_helper'

class DetailedGeoLocationTest < ActiveSupport::TestCase
  test 'cache_key' do
    assert_equal('Oceania/Australia/Northern Territory/Darwin', detailed_geo_locations(:darwin).cache_key)
  end
end
