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
end
