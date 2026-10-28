# frozen_string_literal: true

module VariantPositions
  # The igv.js diagram for one 0-based reference position: a BED track per primer overlapping it, most impacted first,
  # and a BED track per variant allele seen there, most frequent first, and a locus that frames the primers. Each
  # track is coloured by impact (see Disruption). Frequencies are over the sequences collected within the last `days`
  # days. Uses the organism's first taxon, like the organism browser.
  class Track # rubocop:disable Metrics/ClassLength -- the variant and primer tracks each need their own data and naming
    PADDING = 25
    CACHE_TTL = 1.hour

    Allele = Struct.new(:ref_start, :ref_end, :variant_type, :variant, :ref, :sequences, keyword_init: true)
    # also: labels of the other primers at the same coordinates, which the track's description lists
    PrimerTrack = Struct.new(:ref_start, :ref_end, :scheme, :label, :impact, :strand, :also, keyword_init: true)

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
        built = allele_counts.map do |(start, stop, type, variant, ref), sequences|
          Allele.new(ref_start: start, ref_end: stop, variant_type: type, variant:, ref:, sequences:)
        end
        built.sort_by { |a| [-a.sequences, a.ref_start, a.variant_type, a.variant] }
      end
    end

    def allele_name(allele)
      impact = disruption(allele) || 'none'
      "#{allele_label(allele)} (#{percent(allele.sequences)}%), disruption #{impact}"
    end

    # HGVS-like label: A>G for a mismatch, ΔAGT for a deletion of AGT, A>AAAA for an insertion after anchor base A.
    # Falls back to the length or a '?' when the reference base was not recorded.
    def allele_label(allele)
      case allele.variant_type
      when 'X' then "#{allele.ref || '?'}>#{allele.variant}"
      when 'D' then allele.ref.present? ? "Δ#{allele.ref}" : "Δ#{allele.variant.to_i} bp"
      else allele.ref.present? ? "#{allele.ref}>#{allele.ref}#{allele.variant}" : "+#{allele.variant}"
      end
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
        groups = primer_positions.group_by { |position| [position.ref_start, position.ref_end] }
        built = groups.values.map { |group| primer_track(group) }
        built.sort_by { |primer| -(primer.impact || -Float::INFINITY) }
      end
    end

    # one track for primers sharing coordinates: the most recent primer set, then the highest oligo id, names it
    def primer_track(group)
      latest, *others = group.sort_by { |p| [p.oligo.primer_set.created_at, p.oligo.id] }.reverse
      PrimerTrack.new(ref_start: latest.ref_start, ref_end: latest.ref_end, scheme: latest.oligo.primer_set.name,
                      label: latest.oligo.label, impact: group.filter_map { |p| primer_impacts[p.oligo.id] }.max,
                      strand: primer_strand(latest), also: other_primer_labels(others))
    end

    def other_primer_labels(others)
      others.map do |p|
        "#{p.oligo.primer_set.name}: #{p.oligo.label} (impact #{impact_text(primer_impacts[p.oligo.id])})"
      end
    end

    # the track's name: the primer scheme (primer set); the primer and its impact are in the description
    def primer_track_name(primer)
      primer.scheme
    end

    # the track's popup text: the primer and its impact, then any other primers at the same coordinates
    # one line for the primer and its impact, then, when others share its coordinates, a line for each as set: primer
    def primer_description(primer)
      lines = ["#{primer.label}, impact #{impact_text(primer.impact)}"]
      lines += ['Same coordinates:', *primer.also] if primer.also.any?
      lines.join("\n")
    end

    # the variant track's name: the allele; its share of sequences and disruption are in the description
    def allele_track_name(allele)
      allele_label(allele)
    end

    def allele_description(allele)
      "#{percent(allele.sequences)}% of sequences in the window (#{allele.sequences}); " \
        "disruption #{disruption(allele) || 'none'}"
    end

    def primer_name(primer)
      "#{primer.scheme} #{primer.label}, impact #{impact_text(primer.impact)}"
    end

    def primer_label(position)
      "#{position.oligo.primer_set.name} #{position.oligo.label}"
    end

    def impact_text(impact)
      impact ? format('%.1f', impact) : 'none'
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
                 .group(:ref_start, :ref_end, :variant_type, :variant, :ref).count
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

    # the primer's direction along the reference, see Strand
    def primer_strand(position)
      primer_strands.fetch(position.oligo.id, '.')
    end

    def primer_strands
      @primer_strands ||= cached(:primer_strands) do
        primer_positions.to_h do |position|
          [position.oligo.id, Strand.of(position.oligo, position.ref_start, position.ref_end, reference_sequence)]
        end
      end
    end

    def reference_sequence
      @reference_sequence = VariantSite.load_reference_sequence(taxon) unless defined?(@reference_sequence)
      @reference_sequence
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
