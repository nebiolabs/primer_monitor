# frozen_string_literal: true

module VariantPositions
  # A primer's direction along the reference: the strand recorded for it, or one inferred from its sequence when none is
  # recorded (the nightly classifier does not infer it, so the diagrams are the only place it is).
  module Strand
    COMPLEMENT = { 'A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A' }.freeze

    # '+' or '-', or '.' when neither the recorded strand nor the reference settles it
    def self.of(oligo, ref_start, ref_end, reference)
      recorded = oligo.strand
      return recorded if %w[+ -].include?(recorded)

      inferred(oligo, ref_start, ref_end, reference) || '.'
    end

    # '+' when the primer reads along the reference, '-' when its reverse complement does; nil unless one of them
    # matches at least 90% of the bases
    def self.inferred(oligo, ref_start, ref_end, reference)
      sequence = oligo.sequence.to_s.upcase.tr('U', 'T')
      return if reference.nil? || sequence.empty?

      span = reference[ref_start...ref_end]
      return if span.nil? || span.length != sequence.length

      forward = identity(span, sequence)
      reverse = identity(span, reverse_complement(sequence))
      return if [forward, reverse].max < 0.9

      forward >= reverse ? '+' : '-'
    end

    def self.identity(expected, actual)
      expected.chars.zip(actual.chars).count { |a, b| a == b } / expected.length.to_f
    end

    def self.reverse_complement(sequence)
      sequence.reverse.chars.map { |base| COMPLEMENT.fetch(base, 'N') }.join
    end
  end
end
