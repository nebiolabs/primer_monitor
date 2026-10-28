# frozen_string_literal: true

module Lineages
  # Where a lineage's sequences were collected within the window, as map points: each position's count of the
  # lineage's sequences, the sequences from that position in total, and the lineage's share of them.
  class Geography
    def self.call(organism:, lineage:, since:)
      new(organism:, lineage:, since:).call
    end

    def initialize(organism:, lineage:, since:)
      @organism = organism
      @lineage = lineage
      @since = since
    end

    def call
      window = FastaRecord.where(organism_taxon_id: organism.organism_taxa.ids, date_collected: since..)
      lineage_records = window.joins(lineage_call: :lineage).where(lineage_calls: { lineage_id: lineage.id })
      {
        lineage_sequences: lineage_records.count,
        sequences: window.count,
        points: GeoPoints.rows(GeoPoints.sums(lineage_records), GeoPoints.sums(window))
      }
    end

    private

    attr_reader :organism, :lineage, :since
  end
end
