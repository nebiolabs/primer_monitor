# frozen_string_literal: true

require 'test_helper'

class OligoTest < ActiveSupport::TestCase
  test 'the truth' do
    assert_not_nil o = Oligo.new(name: 'test', sequence: 'GCCATATAT')
    assert_equal('test: GCCATATAT', o.to_s)
  end

  test 'fills blank locus and category from the name' do
    o = Oligo.new(name: 'nCoV-2019_7_LEFT_alt0', sequence: 'GCCATATAT', primer_set: primer_sets(:one))

    assert_predicate o, :valid?
    assert_equal %w[7 Forward], [o.locus, o.category]
  end

  test 'keeps locus and category entered by hand' do
    o = Oligo.new(name: 'N1-F', locus: 'X', category: 'Probe', sequence: 'GCCATATAT', primer_set: primer_sets(:one))
    o.valid?

    assert_equal %w[X Probe], [o.locus, o.category]
  end

  test 'label adds locus and category when present' do
    assert_equal 'N1-F (N1 Forward)', Oligo.new(name: 'N1-F').tap(&:valid?).label
    assert_equal 'probe1', Oligo.new(name: 'probe1').label
  end
end
