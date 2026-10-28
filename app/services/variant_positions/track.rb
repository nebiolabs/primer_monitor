# frozen_string_literal: true

module VariantPositions
  # The igv.js diagram for one 0-based reference position: a BED track per variant allele seen there, most frequent
  # first, the primers overlapping the position, and a locus that frames those primers. Frequencies are over the
  # sequences collected within the last `days` days. Uses the organism's first taxon, like the organism browser.
  class Track
    PADDING = 25
    CACHE_TTL = 1.hour

    Allele = Struct.new(:ref_start, :ref_end, :variant_type, :variant, :sequences, keyword_init: true)

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
      "#{allele.variant_type}:#{allele.variant} (#{percent(allele.sequences)}%)"
    end

    # BED line of one allele, by its index in #alleles; nil when there is no such allele. The score is the number of
    # sequences carrying the allele.
    def variant_bed(index)
      allele = alleles[index]
      return if allele.nil?

      bed([[allele.ref_start, [allele.ref_end, allele.ref_start + 1].max, allele_name(allele), allele.sequences]])
    end

    # BED lines of the primers overlapping the position, named by primer set and primer. The score is the number of
    # distinct variant sites in the window overlapping the primer.
    def primer_bed
      counts = primer_variant_counts
      bed(primer_positions.map do |position|
        oligo = position.oligo
        [position.ref_start, position.ref_end, "#{oligo.primer_set.name} #{oligo.label}",
         counts.fetch(oligo.id, 0), strand(oligo)]
      end)
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

    # oligo id => distinct variant sites in the window overlapping that primer, in one query against the nightly
    # overlap view (which is why the page caps its window at the view's 168 days)
    def primer_variant_counts
      @primer_variant_counts ||= cached(:primer_variant_counts) do
        OligoVariantOverlap.where(oligo_id: primer_positions.map { |p| p.oligo.id })
                           .where(date_collected: since..).group(:oligo_id).distinct.count(:variant_id)
      end
    end

    # the diagram makes one request per allele track, all for the same data, so share it for a while
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

    # the primer's direction along the reference; '.' when it is not recorded
    def strand(oligo)
      %w[+ -].include?(oligo.strand) ? oligo.strand : '.'
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
