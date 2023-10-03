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

  test 'python_sam_parser' do
    input_sam = Rails.root.join('test/fixtures/complex_cigar.sam')
    parser = Rails.root.join('lib/parse_alignments.py')
    output = `python #{parser} < #{input_sam}`
    assert_equal "hCoV-19/England/MILK-31C5A2D/2022|EPI_ISL_8755092|2022-01-07\t28362\tD\t---------",
                 output.split("\n")[-1]
  end

  test 'soft clips' do
    input_sam = Rails.root.join('test/fixtures/soft_clip_cigar.sam')
    parser = Rails.root.join('lib/parse_alignments.py')
    output = `python #{parser} < #{input_sam}`
    assert_equal "hCoV-19/USA/CO-CDPHE-2100422558/2021\t29197\tX\tT",
                 output.split("\n")[-1]
  end
end
