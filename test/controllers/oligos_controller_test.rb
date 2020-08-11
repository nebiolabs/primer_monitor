require 'test_helper'

class OligosControllerTest < ActionDispatch::IntegrationTest
  setup do
    @oligo = oligos(:one)
  end

  test "should get index" do
    get oligos_url
    assert_response :success
  end

  test "should get new" do
    get new_oligo_url
    assert_response :success
  end

  test "should create oligo" do
    assert_difference('Oligo.count') do
      post oligos_url, params: { oligo: { name: @oligo.name, sequence: @oligo.sequence } }
    end

    assert_redirected_to oligo_url(Oligo.last)
  end

  test "should show oligo" do
    get oligo_url(@oligo)
    assert_response :success
  end

  test "should get edit" do
    get edit_oligo_url(@oligo)
    assert_response :success
  end

  test "should update oligo" do
    patch oligo_url(@oligo), params: { oligo: { name: @oligo.name, sequence: @oligo.sequence } }
    assert_redirected_to oligo_url(@oligo)
  end

  test "should destroy oligo" do
    assert_difference('Oligo.count', -1) do
      delete oligo_url(@oligo)
    end

    assert_redirected_to oligos_url
  end
end
