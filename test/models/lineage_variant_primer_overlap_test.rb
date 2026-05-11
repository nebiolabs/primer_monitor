# frozen_string_literal: true

require 'test_helper'

class LineageVariantPrimerOverlapTest < ActiveSupport::TestCase
  test 'organism has pct_cutoff' do
    assert_respond_to organisms(:sars_cov2), :pct_cutoff
    assert_equal 1.0, organisms(:sars_cov2).pct_cutoff
  end
end
