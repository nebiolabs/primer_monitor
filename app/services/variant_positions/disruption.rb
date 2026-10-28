# frozen_string_literal: true

module VariantPositions
  # How much one variant would disrupt one primer, ported from lib/visualization/classify_overlaps.py, which the
  # primer-impact tracks are built with. A variant scores each base it covers inside the primer by how far that base is
  # from the primer's 5' and 3' ends, multiplied for indels; the variant's contribution to the primer is 2 ** score.
  # Scores are returned as the exponent (the score).
  module Disruption
    NEAR_5P_SCORE = 1
    NEAR_3P_SCORE = 2
    MID_3P_SCORE = 4
    OTHER_SCORE = 3
    INDEL_MULTIPLIER = 10

    NO_IMPACT_COLOR = '#9e9e9e'
    LOW_RGB = [255, 224, 102].freeze # yellow, for the lowest impact
    HIGH_RGB = [193, 18, 31].freeze # red, for the highest
    SCALE_TOP = 20 # scores at or above this are drawn in the highest colour

    # variant_type is 'X' (mismatch), 'D' (deletion) or 'I' (insertion); spans are 0-based reference ranges with an
    # exclusive end. Returns nil when the variant does not affect the primer.
    def self.primer_score(variant_type:, variant_span:, primer_span:, strand:)
      primer_start = primer_span.begin
      primer_end = primer_span.end
      # an insertion counts only when it starts inside the primer
      return if variant_type == 'I' && !variant_span.begin.between?(primer_start, primer_end - 1)

      overlap = [primer_start, variant_span.begin].max...[primer_end, variant_span.end].min
      return if overlap.size <= 0

      indel = %w[D I].include?(variant_type)
      dist_5p, dist_3p, step = end_distances(overlap.begin, primer_start, primer_end, strand)

      overlap.sum do
        base = base_score(dist_5p, dist_3p)
        dist_5p += step
        dist_3p -= step
        indel ? base * INDEL_MULTIPLIER : base
      end
    end

    # distances of a position from the primer's 5' and 3' ends, and the direction those distances move along the primer
    def self.end_distances(position, primer_start, primer_end, strand)
      if strand == '+'
        [position - primer_start, primer_end - position, 1]
      else
        [primer_end - position, position - primer_start, -1]
      end
    end

    def self.base_score(dist_5p, dist_3p)
      if dist_5p <= 3 then NEAR_5P_SCORE
      elsif dist_3p <= 2 then NEAR_3P_SCORE
      elsif dist_3p.between?(3, 5) then MID_3P_SCORE
      else OTHER_SCORE
      end
    end

    # log2(2 ** e1 + 2 ** e2 + ...), computed without overflowing for large exponents
    def self.log2_sum(exponents)
      top = exponents.max
      top + Math.log2(exponents.sum { |exponent| 2.0**(exponent - top) })
    end

    # track colour for a disruption score: grey when nil (no primer affected), otherwise yellow (low) to red (high)
    def self.color(score)
      return NO_IMPACT_COLOR if score.nil?

      fraction = [score, SCALE_TOP].min / SCALE_TOP.to_f
      red, green, blue = LOW_RGB.zip(HIGH_RGB).map { |low, high| (low + ((high - low) * fraction)).round }
      format('#%<red>02x%<green>02x%<blue>02x', red:, green:, blue:)
    end
  end
end
