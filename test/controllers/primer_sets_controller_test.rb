require 'test_helper'

class PrimerSetsControllerTest < ActionDispatch::IntegrationTest
  setup do
    @primer_set = primer_sets(:one)
  end

  test "should get index" do
    get primer_sets_url
    assert_response :success
  end

  test "should get new" do
    get new_primer_set_url
    assert_response :success
  end

  test "should create primer_set" do
    assert_difference('PrimerSet.count') do
      post primer_sets_url, params: { primer_set: { creator_id: @primer_set.creator_id, name: @primer_set.name } }
    end

    assert_redirected_to primer_set_url(PrimerSet.last)
  end

  test "should show primer_set" do
    get primer_set_url(@primer_set)
    assert_response :success
  end

  test "should get edit" do
    get edit_primer_set_url(@primer_set)
    assert_response :success
  end

  test "should update primer_set" do
    patch primer_set_url(@primer_set), params: { primer_set: { creator_id: @primer_set.creator_id, name: @primer_set.name } }
    assert_redirected_to primer_set_url(@primer_set)
  end

  test "should destroy primer_set" do
    assert_difference('PrimerSet.count', -1) do
      delete primer_set_url(@primer_set)
    end

    assert_redirected_to primer_sets_url
  end
end
