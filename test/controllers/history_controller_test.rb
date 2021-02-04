require "test_helper"

class HistoryControllerTest < ActionDispatch::IntegrationTest
  test "should get show" do
    get history_show_url
    assert_response :success
  end
end
