# frozen_string_literal: true

require 'test_helper'

class RoleTest < ActiveSupport::TestCase
  test 'instantiation' do
    assert_not_nil r = Role.new(name: 'king')
    assert_equal 'king', r.to_s
  end
end
