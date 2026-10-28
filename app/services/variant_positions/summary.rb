# frozen_string_literal: true

module VariantPositions
  # Where sequences collected within the last `days` days carry a variant at one 0-based reference position, and
  # which named lineages those sequences belong to. Fractions are over all sequences in the same location or lineage
  # within the window, the same denominator the notification emails use.
  class Summary
    Row = Struct.new(:name, :variant_count, :sequence_count, :fraction, :latitude, :longitude, keyword_init: true)

    def self.call(organism:, coordinate:, days:)
      new(organism:, coordinate:, days:).call
    end

    def initialize(organism:, coordinate:, days:)
      @organism = organism
      @coordinate = coordinate
      @since = Date.current - days
    end

    def call
      {
        variant_sequences: variant_records.count,
        sequences: window_records.count,
        locations: locations,
        lineages: lineages
      }
    end

    private

    attr_reader :organism, :coordinate, :since

    def taxon_ids
      @taxon_ids ||= organism.organism_taxa.ids
    end

    def window_records
      FastaRecord.where(organism_taxon_id: taxon_ids, date_collected: since..)
    end

    def variant_records
      window_records.where(id: VariantSite.where(organism_taxon_id: taxon_ids).at_position(coordinate)
                                          .select(:fasta_record_id))
    end

    def locations
      totals = window_records.where.not(detailed_geo_location_id: nil).group(:detailed_geo_location_id).count
      variants = variant_records.where.not(detailed_geo_location_id: nil).group(:detailed_geo_location_id).count
      names = DetailedGeoLocation.includes(:detailed_geo_location_alias).where(id: variants.keys).index_by(&:id)

      rows = variants.map do |location_id, variant_count|
        location = names.fetch(location_id)
        alias_record = location.detailed_geo_location_alias
        build_row(location_name(location), variant_count, totals.fetch(location_id, 0),
                  latitude: alias_record&.latitude, longitude: alias_record&.longitude)
      end
      sort_rows(rows)
    end

    def lineages
      totals = named_lineage_counts(window_records)
      variants = named_lineage_counts(variant_records)
      names = Lineage.where(id: variants.keys).index_by(&:id)

      rows = variants.map do |lineage_id, variant_count|
        build_row(names.fetch(lineage_id).name, variant_count, totals.fetch(lineage_id, 0))
      end
      sort_rows(rows)
    end

    def named_lineage_counts(records)
      records.joins(:lineage_call).where.not(lineage_calls: { lineage_id: nil })
             .group('lineage_calls.lineage_id').count
    end

    # matches how the notification emails name a location, without the empty trailing levels
    def location_name(location)
      [location.region, location.subregion, location.division, location.subdivision].compact.join('/')
    end

    def sort_rows(rows)
      rows.sort_by { |row| [-row.variant_count, row.name] }
    end

    def build_row(name, variant_count, sequence_count, latitude: nil, longitude: nil)
      Row.new(name:, variant_count:, sequence_count:,
              fraction: sequence_count.zero? ? 0.0 : variant_count.to_f / sequence_count,
              latitude:, longitude:)
    end
  end
end
