# frozen_string_literal: true

require 'test_helper'

class WelcomeControllerTest < ActionDispatch::IntegrationTest
  test 'should get landing page' do
    get '/'
    assert_response :success
  end
end
