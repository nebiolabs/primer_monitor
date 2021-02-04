# frozen_string_literal: true

require 'test_helper'

class WelcomeControllerTest < ActionDispatch::IntegrationTest
  test 'should get show' do
    get welcome_index_url
    assert_response :success
  end
end
