# frozen_string_literal: true

require 'test_helper'

class OrganismsControllerTest < ActionDispatch::IntegrationTest
  setup do
    @organism = organisms(:sars_cov2)
    sign_in(users(:admin_user))
  end

  test 'should get show' do
    get organisms_url
    assert_response :success
  end

  test 'should get new' do
    get new_organism_url
    assert_response :success
  end

  test 'should create organism' do
    assert_difference('Organism.count') do
      post organisms_url,
           params: { organism: { name: 'test organism', slug: 'test',
                                 organism_taxa: [organism_taxa(:sars_cov_2_taxon)] } }
    end

    assert_redirected_to organism_url(Organism.last)
  end

  test 'should show organism' do
    Organism.any_instance.stubs(:primer_sets_config).returns([{}, {}])
    get organism_url(@organism)
    assert_response :success
  end

  test 'should get edit' do
    get edit_organism_url(@organism)
    assert_response :success
  end

  test 'should update organism' do
    patch organism_url(@organism), params:
      { organism: { name: '${@organism.name}∆' } }
    assert_redirected_to organism_url(@organism)
  end

  test 'should destroy organism' do
    assert_difference('Organism.count', -1) do
      delete organism_url(@organism)
    end

    assert_redirected_to organisms_url
  end
end
