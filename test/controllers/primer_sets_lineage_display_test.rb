# frozen_string_literal: true

require 'test_helper'

# the primer set page shows the primers against the variants of one lineage, by default the most prevalent recent one
class PrimerSetsLineageDisplayTest < ActionDispatch::IntegrationTest
  setup do
    @primer_set = primer_sets(:one)
    sign_in(users(:admin_user))
  end

  test 'show offers the lineages and selects the most prevalent recent one' do
    stub_lineage_data(frequencies: { 'B' => 3, 'XBB' => 12, 'all' => 15 })
    get primer_set_url(@primer_set)

    assert_select '#primer_set_lineage_select option', count: 3
    assert_select '#primer_set_lineage_select option[selected]', count: 1 do |option|
      assert_equal 'XBB', option.first['value']
    end
    assert_select '#variant_table_section'
  end

  test 'show selects the lineage given in the URL' do
    stub_lineage_data(frequencies: { 'B' => 3, 'XBB' => 12, 'all' => 15 })
    get primer_set_url(@primer_set, lineage: 'B')

    assert_select '#primer_set_lineage_select option[selected]' do |option|
      assert_equal 'B', option.first['value']
    end
  end

  test 'show leaves out the lineage picker when the lineage data is unavailable' do
    stub_lineage_data(fetched: false)
    get primer_set_url(@primer_set)

    assert_select '#igv'
    assert_select '#primer_set_lineage_select', count: 0
  end

  private

  def stub_lineage_data(fetched: true, frequencies: {})
    @primer_set.update!(status: :complete)
    Organism.any_instance.stubs(:primer_sets_config)
            .returns([{ data_server: 'http://data', organism_slug: 'sars-cov-2' }, { @primer_set.name => 'Charite' }])
    track_data = if fetched
                   { data_fetched: true, primer_sets: { 'Charite' => @primer_set.name }, default_tracks: ['Charite'],
                     lineage_sets: { 'B' => 'B*', 'XBB' => 'XBB*', 'all' => 'All' }, default_lineage: 'all' }
                 else
                   { data_fetched: false }
                 end
    Organism.any_instance.stubs(:lineage_variants_data).returns(track_data)
    Organism.any_instance.stubs(:grouped_lineage_frequencies).returns(frequencies)
  end
end
