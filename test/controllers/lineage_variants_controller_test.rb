# frozen_string_literal: true

require 'test_helper'

class LineageVariantsControllerTest < ActionDispatch::IntegrationTest
  setup do
    @organism = organisms(:sars_cov2)
    sign_in(users(:admin_user))
  end

  test 'variant_overlaps returns 400 without lineage param' do
    get variant_overlaps_organism_lineage_variants_url(
      organism_slug: @organism.slug,
      primer_sets: ['Charité'],
      format: :json
    )
    assert_response :bad_request
  end

  test 'variant_overlaps returns 400 without primer_sets param' do
    get variant_overlaps_organism_lineage_variants_url(
      organism_slug: @organism.slug,
      lineage: 'XBB',
      format: :json
    )
    assert_response :bad_request
  end

  test 'variant_overlaps returns 404 for unknown organism' do
    get variant_overlaps_organism_lineage_variants_url(
      organism_slug: 'no-such-organism',
      lineage: 'XBB',
      primer_sets: ['Charité'],
      format: :json
    )
    assert_response :not_found
  end

  test 'variant_overlaps returns JSON for valid params' do
    get variant_overlaps_organism_lineage_variants_url(
      organism_slug: @organism.slug,
      lineage: 'XBB',
      primer_sets: ['Charité'],
      format: :json
    )
    assert_response :success
    data = JSON.parse(response.body)
    assert data.key?('variants')
  end

  test 'variant_overlaps returns matching variants with oligos' do
    get variant_overlaps_organism_lineage_variants_url(
      organism_slug: @organism.slug,
      lineage: 'XBB',
      primer_sets: ['Charité'],
      format: :json
    )
    assert_response :success
    data = JSON.parse(response.body)
    assert_equal 1, data['variants'].length

    v = data['variants'].first
    assert_equal 100, v['ref_start']
    assert_equal 101, v['ref_end']
    assert_equal 'X', v['variant_type']
    assert_equal 'T', v['variant']
    assert_in_delta 15.5, v['frequency_pct'], 0.01

    assert_equal 1, v['oligos'].length
    o = v['oligos'].first
    assert_equal 'probe1', o['name']
    assert_equal 'GCAATTTATATACATATA', o['sequence']
    assert_equal 'Charité', o['primer_set']
  end

  test 'variant_overlaps excludes variants for other primer sets' do
    get variant_overlaps_organism_lineage_variants_url(
      organism_slug: @organism.slug,
      lineage: 'XBB',
      primer_sets: ['CDC'],
      format: :json
    )
    assert_response :success
    data = JSON.parse(response.body)
    assert_equal 0, data['variants'].length
  end

  test 'variant_overlaps excludes variants for other lineage groups' do
    get variant_overlaps_organism_lineage_variants_url(
      organism_slug: @organism.slug,
      lineage: 'JN.1',
      primer_sets: ['Charité'],
      format: :json
    )
    assert_response :success
    data = JSON.parse(response.body)
    assert_equal 0, data['variants'].length
  end

  test 'variant_overlaps includes first_seen and last_seen keys' do
    get variant_overlaps_organism_lineage_variants_url(
      organism_slug: @organism.slug,
      lineage: 'XBB',
      primer_sets: ['Charité'],
      format: :json
    )
    assert_response :success
    v = JSON.parse(response.body)['variants'].first
    assert v.key?('first_seen'), 'missing first_seen'
    assert v.key?('last_seen'),  'missing last_seen'
    assert v['first_seen'].key?('date')
    assert v['first_seen'].key?('lineage')
    assert v['first_seen'].key?('location')
  end
end
