# frozen_string_literal: true

require 'test_helper'

class UsersControllerTest < ActionDispatch::IntegrationTest
  setup do
    @user = users(:admin_user)
    sign_in(@user)
  end

  test 'should get index' do
    get users_url
    assert_response :success
  end

  test 'should get new' do
    sign_out(@user)
    get new_user_registration_url
    assert_response :success
  end

  test 'should create user' do
    sign_out(@user)
    assert_difference('User.count') do
      post user_registration_url, params: { user: { email: 'test@example.org', first: 'test', last: 'user',
                                                    password: 'test123', password_confirmation: 'test123' } }
    end

    assert_redirected_to root_url
  end

  test 'should show user' do
    get user_url(@user)
    assert_response :success
  end

  test 'should get edit' do
    get edit_user_registration_path(@user)
    assert_response :success
  end

  test 'should update user' do
    patch user_url(@user), params: { user: { email: @user.email, first: @user.first, last: @user.last } }
    assert_redirected_to user_url(@user)
  end

  test 'should destroy user' do
    assert_difference('User.count', -1) do
      delete user_url(@user)
    end

    assert_redirected_to users_url
  end
end
