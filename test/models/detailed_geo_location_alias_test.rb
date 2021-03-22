# frozen_string_literal: true

require 'test_helper'

class DetailedGeoLocationAliasTest < ActiveSupport::TestCase
  test 'instantiation' do
    assert_not_nil DetailedGeoLocationAlias.new
  end

  test 'subscribable_geolocations' do
    assert_empty DetailedGeoLocationAlias.subscribable
    orig_rec = FastaRecord.first
    20.times.each do |i|
      rec = orig_rec.dup
      rec.strain += i.to_s
      rec.save!
    end
    assert_includes(DetailedGeoLocationAlias.subscribable, orig_rec.detailed_geo_location.detailed_geo_location_alias)
  end
end
