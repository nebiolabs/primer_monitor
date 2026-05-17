# frozen_string_literal: true

module LineageVariants
  class OverlapResponseBuilder
    def self.call(organism:, lineage_group_key:, primer_set_names:)
      new(organism:, lineage_group_key:, primer_set_names:).variants
    end

    def initialize(organism:, lineage_group_key:, primer_set_names:)
      @organism = organism
      @lineage_group_key = lineage_group_key
      @primer_set_names = primer_set_names
    end

    def variants
      grouped_variants.values
    end

    private

    attr_reader :organism, :lineage_group_key, :primer_set_names

    def grouped_variants
      fetch_rows.each_with_object({}) do |row, variants_map|
        key = [row.ref_start, row.ref_end, row.variant_type, row.variant]
        variants_map[key] ||= grouped_variant(row)
        variants_map[key][:oligos] << grouped_oligo(row)
      end
    end

    def fetch_rows
      LineageVariantPrimerOverlap.for_variant_overlaps(
        organism_id: organism.id,
        lineage_group_key:,
        primer_set_names:
      )
    end

    def grouped_variant(row)
      {
        ref_start: row.ref_start,
        ref_end: row.ref_end,
        variant_type: row.variant_type,
        variant: row.variant,
        ref: row.ref,
        frequency_pct: row.frequency_pct,
        first_seen: { date: row.first_seen_date, lineage: row.first_seen_lineage, location: row.first_seen_location },
        last_seen: { date: row.last_seen_date, lineage: row.last_seen_lineage, location: row.last_seen_location },
        oligos: []
      }
    end

    def grouped_oligo(row)
      {
        name: row.oligo_name,
        sequence: row.oligo_sequence,
        oligo_start: row.oligo_start,
        oligo_end: row.oligo_end,
        strand: row.oligo_strand,
        primer_set: row.primer_set_name
      }
    end
  end
end
