# frozen_string_literal: true

require 'test_helper'

class FastaRecordTest < ActiveSupport::TestCase
  test 'metadata_parse' do
    new_records = FastaRecord.parse(Rails.root.join('test/fixtures/metadata.tsv'))
    assert_equal(7, new_records.size) # AUS/NT01/2020 and AUS/NT02/2020 are already in via the fixtures
    assert(new_records[0].valid?, new_records[0].errors.full_messages)
    assert(new_records[0].detailed_geo_location.valid?, new_records[0].detailed_geo_location.errors.full_messages)
    assert(new_records[0].detailed_geo_location
                         .detailed_geo_location_alias.valid?,
           new_records[0].detailed_geo_location
                         .detailed_geo_location_alias.errors.full_messages)
  end
end
