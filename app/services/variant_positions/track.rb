# frozen_string_literal: true

module VariantPositions
  # The igv.js diagram for one 0-based reference position: a BED track per primer overlapping it, most impacted first,
  # and a BED track per variant allele seen there, most frequent first, and a locus that frames the primers. Each
  # track is coloured by impact (see Disruption). Frequencies are over the sequences collected within the last `days`
  # days. Uses the organism's first taxon, like the organism browser.
  class Track # rubocop:disable Metrics/ClassLength -- the variant and primer tracks each need their own data and naming
    PADDING = 25
    CACHE_TTL = 1.hour

    Allele = Struct.new(:ref_start, :ref_end, :variant_type, :variant, :sequences, keyword_init: true)
    COMPLEMENT = { 'A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A' }.freeze
    PrimerTrack = Struct.new(:ref_start, :ref_end, :name, :impact, :strand, keyword_init: true)

    def initialize(organism:, coordinate:, days:)
      @coordinate = coordinate
      @since = Date.current - days
      @taxon = organism.organism_taxa.first
    end

    def reference_accession
      taxon&.reference_accession
    end

    # the variant alleles at the position, most frequent first; each one is a track of its own so the order holds
    def alleles
      @alleles ||= cached(:alleles) do
        built = allele_counts.map do |(start, stop, type, variant), sequences|
          Allele.new(ref_start: start, ref_end: stop, variant_type: type, variant:, sequences:)
        end
        built.sort_by { |a| [-a.sequences, a.ref_start, a.variant_type, a.variant] }
      end
    end

    def allele_name(allele)
      impact = disruption(allele) || 'none'
      "#{allele.variant_type}:#{allele.variant} (#{percent(allele.sequences)}%), disruption #{impact}"
    end

    # the largest disruption the allele causes to a primer overlapping the position (see Disruption), or nil when
    # it affects none of them
    def disruption(allele)
      variant_span = allele.ref_start...[allele.ref_end, allele.ref_start + 1].max
      primer_positions.filter_map do |position|
        Disruption.primer_score(variant_type: allele.variant_type, variant_span:,
                                primer_span: position.ref_start...position.ref_end, strand: primer_strand(position))
      end.max
    end

    # track colour for the allele, see Disruption.color
    def allele_color(allele)
      Disruption.color(disruption(allele))
    end

    # BED line of one allele, by its index in #alleles; nil when there is no such allele. The score is the number of
    # sequences carrying the allele.
    def variant_bed(index)
      allele = alleles[index]
      return if allele.nil?

      bed([[allele.ref_start, [allele.ref_end, allele.ref_start + 1].max, allele_name(allele), allele.sequences]])
    end

    # the primers overlapping the position, most impacted first; each is a track of its own so its colour holds
    def primers
      @primers ||= begin
        built = primer_positions.map do |position|
          oligo = position.oligo
          PrimerTrack.new(ref_start: position.ref_start, ref_end: position.ref_end,
                          name: "#{oligo.primer_set.name} #{oligo.label}", impact: primer_impacts[oligo.id],
                          strand: primer_strand(position))
        end
        built.sort_by { |primer| -(primer.impact || -Float::INFINITY) }
      end
    end

    def primer_name(primer)
      "#{primer.name}, impact #{primer.impact ? format('%.1f', primer.impact) : 'none'}"
    end

    # track colour for the primer, from its impact, see Disruption.color
    def primer_color(primer)
      Disruption.color(primer.impact)
    end

    # BED line of one primer, by its index in #primers; nil when there is no such primer. The score is the rounded
    # impact.
    def primer_bed(index)
      primer = primers[index]
      return if primer.nil?

      bed([[primer.ref_start, primer.ref_end, primer_name(primer), (primer.impact || 0).round, primer.strand]])
    end

    # 1-based igv locus framing the overlapping primers and the position
    def locus
      spans = primer_positions.map { |p| [p.ref_start, p.ref_end] }
      spans = [[coordinate, coordinate + 1]] if spans.empty?
      start = [spans.map(&:first).min - PADDING, 0].max
      stop = spans.map(&:last).max + PADDING
      "#{reference_accession}:#{start + 1}-#{stop}"
    end

    private

    attr_reader :coordinate, :since, :taxon

    def allele_counts
      VariantSite.where(organism_taxon_id: taxon&.id).at_position(coordinate)
                 .group(:ref_start, :ref_end, :variant_type, :variant).count
    end

    # share of the sequences in the window that carry the allele
    def percent(count)
      total = window_sequence_count
      total.zero? ? 0.0 : (count * 100.0 / total).round(1)
    end

    def window_sequence_count
      @window_sequence_count ||= cached(:sequence_count) do
        taxon ? FastaRecord.where(organism_taxon_id: taxon.id, date_collected: since..).count : 0
      end
    end

    # oligo id => the primer's impact: log2 of the sum over its distinct variant sites in the window of 2 ** disruption,
    # the composite score the classifier thresholds. Nil for a primer no variant in the window affects. Read from the
    # nightly overlap view, which is why the page caps its window at the view's 168 days.
    def primer_impacts
      @primer_impacts ||= cached(:primer_impacts) do
        spans = primer_positions.index_by { |position| position.oligo.id }
        scores = Hash.new { |hash, key| hash[key] = [] }
        overlap_rows(spans.keys).each do |oligo_id, type, start, stop|
          score = primer_disruption(spans.fetch(oligo_id), type, start...stop)
          scores[oligo_id] << score if score
        end
        scores.transform_values { |exponents| Disruption.log2_sum(exponents) }
      end
    end

    def overlap_rows(oligo_ids)
      OligoVariantOverlap.where(oligo_id: oligo_ids).where(date_collected: since..).distinct
                         .pluck(:oligo_id, :variant_type, :variant_start, :variant_end)
    end

    def primer_disruption(position, variant_type, variant_span)
      Disruption.primer_score(variant_type:, variant_span:, primer_span: position.ref_start...position.ref_end,
                              strand: primer_strand(position))
    end

    # the diagram makes one request per track, all for the same data, so share it for a while
    def cached(name, &)
      Rails.cache.fetch(['variant_position_track', taxon&.id, coordinate, since, name], expires_in: CACHE_TTL, &)
    end

    def primer_positions
      return [] if taxon.nil?

      @primer_positions ||= OligoAlignmentPosition.where(organism_taxon: taxon)
                                                  .where('ref_start <= :p AND ref_end > :p', p: coordinate)
                                                  .includes(oligo: :primer_set)
                                                  .order(:ref_start)
                                                  .to_a
    end

    # the primer's direction along the reference: the recorded strand, or one inferred from its sequence when none is
    # recorded; '.' when neither settles it
    def primer_strand(position)
      primer_strands.fetch(position.oligo.id, '.')
    end

    def primer_strands
      @primer_strands ||= cached(:primer_strands) do
        primer_positions.to_h { |position| [position.oligo.id, strand_of(position)] }
      end
    end

    def strand_of(position)
      recorded = position.oligo.strand
      return recorded if %w[+ -].include?(recorded)

      inferred_strand(position) || '.'
    end

    # '+' when the primer reads along the reference, '-' when its reverse complement does; nil unless one of them
    # matches at least 90% of the bases
    def inferred_strand(position)
      reference = reference_sequence
      sequence = position.oligo.sequence.to_s.upcase.tr('U', 'T')
      return if reference.nil? || sequence.empty?

      span = reference[position.ref_start...position.ref_end]
      return if span.nil? || span.length != sequence.length

      forward = identity(span, sequence)
      reverse = identity(span, reverse_complement(sequence))
      return if [forward, reverse].max < 0.9

      forward >= reverse ? '+' : '-'
    end

    def reference_sequence
      @reference_sequence = VariantSite.load_reference_sequence(taxon) unless defined?(@reference_sequence)
      @reference_sequence
    end

    def identity(expected, actual)
      expected.chars.zip(actual.chars).count { |a, b| a == b } / expected.length.to_f
    end

    def reverse_complement(sequence)
      sequence.reverse.chars.map { |base| COMPLEMENT.fetch(base, 'N') }.join
    end

    # igv.js splits BED lines on whitespace unless told otherwise, so names must not contain any
    def bed(lines)
      rows = lines.map do |start, stop, name, score, strand = '.'|
        [reference_accession, start, stop, name.gsub(/\s+/, '_'), score, strand].join("\t")
      end
      rows.join("\n")
    end
  end
end
