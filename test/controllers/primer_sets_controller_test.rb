# frozen_string_literal: true

require 'test_helper'

class PrimerSetsControllerTest < ActionDispatch::IntegrationTest
  setup do
    @primer_set = primer_sets(:one)
    sign_in(users(:admin_user))
  end

  test 'should get index' do
    Organism.any_instance.stubs(:primer_sets_config).returns([{}, {}])
    get organism_primer_sets_url(organisms(:sars_cov2))
    assert_response :success
  end

  test 'should get new' do
    get new_primer_set_url
    assert_response :success
  end

  test 'should create primer_set' do
    assert_difference('PrimerSet.count') do
      post primer_sets_url, params: { primer_set: {
        user_id: @primer_set.user_id,
        name: 'test',
        organism_id: organisms(:sars_cov2).id,
        amplification_method_id: amplification_methods(:qPCR).id,
        oligos_attributes: [{ name: 'F', sequence: 'ATCG' },
                            { name: 'R', sequence: 'CGTA' }]
      } }
    end

    assert_redirected_to edit_primer_set_url(PrimerSet.last)
  end

  test 'should show primer_set' do
    Organism.any_instance.stubs(:primer_sets_config).returns([{}, {}])
    get primer_set_url(@primer_set)
    assert_response :success
  end

  test 'should get edit' do
    get edit_primer_set_url(@primer_set)
    assert_response :success
  end

  test 'should update primer_set' do
    patch primer_set_url(@primer_set), params:
      { primer_set: { user_id: @primer_set.user_id, name: "#{@primer_set.name}∆" } }
    assert_redirected_to edit_primer_set_url(@primer_set)
  end

  test 'should destroy primer_set' do
    assert_difference('PrimerSet.count', -1) do
      delete primer_set_url(@primer_set)
    end

    assert_redirected_to primer_sets_url
  end
end
