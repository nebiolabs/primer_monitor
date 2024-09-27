# frozen_string_literal: true

require 'test_helper'

class OrganismsControllerTest < ActionDispatch::IntegrationTest
  include Devise::Test::IntegrationHelpers

  setup do
    sign_in(users(:admin_user))
    @organism = organisms(:example)
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
      post organisms_url, params: { organism: { name: @organism.name } }
    end

    assert_redirected_to organism_url(Organism.last)
  end

  test 'should show organism' do
    get organism_url(@organism.slug)
    assert_response :success
  end

  test 'should get edit' do
    get edit_organism_url(@organism.slug)
    assert_response :success
  end

  test 'should update organism' do
    patch organism_url(@organism.slug), params:
      { organism: { name: @organism.name } }
    assert_redirected_to organism_url(@organism.slug)
  end

  test 'should destroy organism' do
    assert_difference('Organism.count', -1) do
      delete organism_url(@organism.slug)
    end

    assert_redirected_to organisms_url
  end
end
