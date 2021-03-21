# frozen_string_literal: true

require 'test_helper'

class OrganismTest < ActiveSupport::TestCase
  test 'instantiation' do
    assert_not_nil Organism.new
  end
end
