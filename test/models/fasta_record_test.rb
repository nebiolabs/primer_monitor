# frozen_string_literal: true

require 'test_helper'

class FastaRecordTest < ActiveSupport::TestCase
  test 'metadata_parse' do
    new_records = FastaRecord.parse(Rails.root.join('test/fixtures/metadata.tsv'))
    assert_equal(7, new_records.size) # AUS/NT01/2020 and AUS/NT02/2020 are already in via the fixtures
  end
end
