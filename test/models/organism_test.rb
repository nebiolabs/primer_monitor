# frozen_string_literal: true

require 'test_helper'

class OrganismTest < ActiveSupport::TestCase
  setup do
    @organism = Organism.new
  end

  test 'instantiation' do
    assert_not_nil @organism
  end

  test 'process_lineage_variants_data returns correct data structure on valid JSON' do
    tracks_json = '{"track1": "value1"}'
    defaults_json = '{"tracks": "default_tracks", "lineage": "default_lineage"}'
    lineages_json = '{"lineage1": "value1"}'

    expected = {
      data_fetched: true,
      primer_sets: {'track1' => 'value1'},
      default_tracks: 'default_tracks',
      lineage_sets: {'lineage1' => 'value1'},
      default_lineage: 'default_lineage'
    }

    assert_equal expected, @organism.process_lineage_variants_data(tracks_json, defaults_json, lineages_json)
  end

  test 'process_lineage_variants_data handles JSON parsing errors gracefully' do
    invalid_json = 'invalid_json'

    assert_raises(JSON::ParserError) do
      @organism.process_lineage_variants_data(invalid_json, invalid_json, invalid_json)
    end
  end

  test 'process_lineage_variants_data returns empty hashes for empty JSON strings' do
    empty_json = '{}'

    expected = {
      data_fetched: true,
      primer_sets: {},
      default_tracks: nil,
      lineage_sets: {},
      default_lineage: nil
    }

    assert_equal expected, @organism.process_lineage_variants_data(empty_json, empty_json, empty_json)
  end
end
