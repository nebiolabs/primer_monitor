# frozen_string_literal: true

require 'test_helper'

class VariantPositionsControllerTest < ActionDispatch::IntegrationTest
  setup do
    @organism = organisms(:sars_cov2)
    taxon = organism_taxa(:sars_cov_2_taxon)
    record = FastaRecord.create!(strain: 'position/1', detailed_geo_location: detailed_geo_locations(:darwin),
                                 date_collected: Date.current, organism_taxon_id: taxon.id)
    VariantSite.create!(fasta_record_id: record.id, organism_taxon_id: taxon.id, ref_start: 100, ref_end: 101,
                        variant_type: 'X', variant: 'A', ref: nil)
    %w[variant_overlaps oligo_variant_overlaps].each do |view|
      ActiveRecord::Base.connection.execute("REFRESH MATERIALIZED VIEW #{view}")
    end
  end

  test 'shows the geographic distribution of the position' do
    get organism_variant_position_url(organism_slug: @organism.slug, coordinate: 100)

    assert_response :success
    assert_select 'td', text: 'Oceania/Australia/Northern Territory/Darwin'
  end

  test 'accepts a window in days' do
    get organism_variant_position_url(organism_slug: @organism.slug, coordinate: 100, days: 7)

    assert_response :success
    assert_select 'h1', text: /position 100/
  end

  test 'returns 404 for an unknown organism' do
    get organism_variant_position_url(organism_slug: 'no-such-organism', coordinate: 100)

    assert_response :not_found
  end

  test 'page carries the diagram config and the map data' do
    get organism_variant_position_url(organism_slug: @organism.slug, coordinate: 100)

    assert_select 'div#variant_igv'
    assert_select 'p#variant_position_config', text: /primer_tracks/
    assert_select 'div#variant_map[data-locations*="Australia"]'
    assert_select 'tr.is-clickable[data-map-point="Australia"]'
  end

  test 'serves the variant and primer tracks as BED' do
    get organism_variant_position_url(organism_slug: @organism.slug, coordinate: 100, track: 'variant', allele: 0,
                                      format: :text)

    assert_response :success
    assert_equal 'text/plain', response.media_type
    assert_includes response.body, "100\t101\tX:A_(100.0%),_disruption_4\t1\t."

    get organism_variant_position_url(organism_slug: @organism.slug, coordinate: 100, track: 'variant', allele: 9,
                                      format: :text)

    assert_response :not_found

    get organism_variant_position_url(organism_slug: @organism.slug, coordinate: 100, track: 'primer', primer: 0,
                                      format: :text)

    assert_includes response.body, '95'
  end
end
