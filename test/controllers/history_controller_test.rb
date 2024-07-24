# frozen_string_literal: true

require 'test_helper'

class HistoryControllerTest < ActionDispatch::IntegrationTest
  test 'should get show' do
    get history_url
    assert_response :success
  end
end
