# frozen_string_literal: true

require 'test_helper'

class LineageVariantPrimerOverlapTest < ActiveSupport::TestCase
  test 'organism has pct_cutoff' do
    assert_respond_to organisms(:sars_cov2), :pct_cutoff
    assert_equal 1.0, organisms(:sars_cov2).pct_cutoff
  end

  test 'valid overlap can be created' do
    overlap = lineage_variant_primer_overlaps(:one)
    assert overlap.valid?
    assert_equal 'XBB', overlap.lineage_group_key
    assert_equal 100, overlap.ref_start
    assert_equal 'X', overlap.variant_type
  end
end
