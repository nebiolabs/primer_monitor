# frozen_string_literal: true

require 'test_helper'

class UserEmailConfirmationsControllerTest < ActionDispatch::IntegrationTest
  test 'should get show' do
    @user = User.first
    get user_email_confirmation_path(@user.perishable_token)
    assert_response :success
  end
end
