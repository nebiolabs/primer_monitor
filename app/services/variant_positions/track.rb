# frozen_string_literal: true

module VariantPositions
  # The igv.js diagram for one 0-based reference position: BED tracks for the variants seen there and for the primers
  # covering it, and a locus that frames those primers. Uses the organism's first taxon, like the organism browser.
  class Track
    PADDING = 25
    FALLBACK_SPAN = 100

    def initialize(organism:, coordinate:)
      @organism = organism
      @coordinate = coordinate
      @taxon = organism.organism_taxa.first
    end

    def reference_accession
      taxon&.reference_accession
    end

    # BED lines of the variants at the position, with how many sequences carry each one
    def variant_bed
      counts = VariantSite.where(organism_taxon_id: taxon&.id).at_position(coordinate)
                          .group(:ref_start, :ref_end, :variant_type, :variant).count
      lines = counts.map do |(start, stop, type, variant), count|
        [start, [stop, start + 1].max, "#{type}:#{variant} (#{count} seqs)"]
      end
      bed(lines)
    end

    # BED lines of the primers covering the position, named by primer set and primer
    def primer_bed
      bed(primer_positions.map do |position|
        [position.ref_start, position.ref_end, "#{position.oligo.primer_set.name} #{position.oligo.label}"]
      end)
    end

    # 1-based igv locus framing the covering primers and the position
    def locus
      spans = primer_positions.map { |p| [p.ref_start, p.ref_end] }
      spans = [[coordinate, coordinate + 1]] if spans.empty?
      start = [spans.map(&:first).min - PADDING, 0].max
      stop = spans.map(&:last).max + PADDING
      "#{reference_accession}:#{start + 1}-#{stop}"
    end

    private

    attr_reader :coordinate, :taxon

    def primer_positions
      return [] if taxon.nil?

      @primer_positions ||= OligoAlignmentPosition.where(organism_taxon: taxon)
                                                  .where('ref_start <= :p AND ref_end > :p', p: coordinate)
                                                  .includes(oligo: :primer_set)
                                                  .order(:ref_start)
                                                  .to_a
    end

    def bed(lines)
      lines.map { |start, stop, name| [reference_accession, start, stop, name, 0, '+'].join("\t") }.join("\n")
    end
  end
end
