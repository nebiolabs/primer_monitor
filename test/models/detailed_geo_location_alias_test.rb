# frozen_string_literal: true

require 'test_helper'

class DetailedGeoLocationAliasTest < ActiveSupport::TestCase
  test 'instantiation' do
    assert_not_nil DetailedGeoLocationAlias.new
  end

  test 'World is always subscribable; a specific location once it has enough sequences' do
    orig_rec = FastaRecord.first
    location = orig_rec.detailed_geo_location.detailed_geo_location_alias

    assert_equal [detailed_geo_location_aliases(:world)], DetailedGeoLocationAlias.subscribable

    20.times.each do |i|
      rec = orig_rec.dup
      rec.strain += i.to_s
      rec.save!
    end

    assert_includes DetailedGeoLocationAlias.subscribable, location
  end
end
