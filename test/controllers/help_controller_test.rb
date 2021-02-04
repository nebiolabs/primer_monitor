require "test_helper"

class HelpControllerTest < ActionDispatch::IntegrationTest
  test "should get show" do
    get help_show_url
    assert_response :success
  end
end
