# frozen_string_literal: true

require 'test_helper'

module VariantPositions
  class DisruptionTest < ActiveSupport::TestCase
    # a primer spanning 95...115 on the forward strand
    def score(type, start, stop, strand: '+')
      Disruption.primer_score(variant_type: type, variant_span: start...stop, primer_span: 95...115, strand:)
    end

    test 'a mismatch scores by its distance from the primer ends' do
      assert_equal 4, score('X', 110, 111) # 3 to 5 bases from the 3' end
      assert_equal 1, score('X', 97, 98)   # within 3 bases of the 5' end
      assert_equal 2, score('X', 113, 114) # within 2 bases of the 3' end
    end

    test 'an indel multiplies each base by 10' do
      assert_equal 120, score('D', 110, 113) # three bases, each scoring 4 * 10
    end

    test 'the reverse strand measures distance from the other end' do
      assert_equal 3, score('X', 110, 111, strand: '-')
      assert_equal 4, score('X', 100, 101, strand: '-')
    end

    test 'variants outside the primer, and insertions that do not start in it, do not affect it' do
      assert_nil score('X', 200, 201)
      assert_nil score('I', 200, 204)
    end
  end
end
