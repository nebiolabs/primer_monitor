# frozen_string_literal: true

require 'test_helper'

class UsersControllerTest < ActionDispatch::IntegrationTest
  include Devise::Test::IntegrationHelpers

  setup do
    # @admin = users(:admin_user)
    @user = users(:one)
  end

  test 'should get show' do
    sign_in @user
    get user_url @user
    assert_response :success
  end

  test 'should get new' do
    get new_user_registration_url
    assert_response :success
  end

  test 'should create user' do
    assert_difference('User.count') do
      post user_registration_url, params: { user: { email: 'test@example.org', first: 'test', last: 'user',
                                                    password: 'test123', password_confirmation: 'test123' } }
    end

    assert_redirected_to root_url
  end

  test 'should show user' do
    sign_in @user
    get user_url(@user)
    assert_response :success
  end

  test 'should get edit' do
    sign_in @user
    get edit_user_registration_path(@user)
    assert_response :success
  end

  test 'should update user' do
    sign_in @user
    patch user_url(@user), params: { user: { email: @user.email, first: @user.first, last: @user.last } }
    assert_redirected_to user_url(@user)
  end

  # TODO: new deactivate user functionality

  #  test 'should destroy user' do
  #  sign_in @user
  #  assert_difference('User.count', -1) do
  #    delete user_registration_url(@user)
  #  end
  #
  #  assert_redirected_to users_url
  # end
end
