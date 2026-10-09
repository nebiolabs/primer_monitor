# frozen_string_literal: true

require 'test_helper'

class HealthCheckTest < ActionDispatch::IntegrationTest
  test '/up reports the app is up without signing in' do
    get rails_health_check_path

    assert_response :success
  end
end
