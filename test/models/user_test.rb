# frozen_string_literal: true

require 'test_helper'

class UserTest < ActiveSupport::TestCase
  test 'instantiation' do
    assert_not_nil User.new
  end
end
