# frozen_string_literal: true

require 'test_helper'

class UserRoleTest < ActiveSupport::TestCase
  test 'user_role instantiation' do
    assert_not_nil UserRole.new
  end
end
