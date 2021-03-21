# frozen_string_literal: true

require 'test_helper'

class OligoTest < ActiveSupport::TestCase
  test 'the truth' do
    assert_not_nil o = Oligo.new(name: 'test', sequence: 'GCCATATAT')
    assert_equal('test: GCCATATAT', o.to_s)
  end
end
