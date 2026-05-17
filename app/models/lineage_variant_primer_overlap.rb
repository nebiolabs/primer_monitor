# frozen_string_literal: true

class LineageVariantPrimerOverlap < ApplicationRecord
  belongs_to :organism
  belongs_to :oligo

  def self.for_variant_overlaps(organism_id:, lineage_group_key:, primer_set_names:)
    joins(oligo: [:primer_set, { oligo_alignment_positions: :organism_taxon }])
      .where(
        lineage_variant_primer_overlaps: {
          organism_id:, lineage_group_key:
        },
        primer_sets: { name: primer_set_names }, organism_taxa: { organism_id: }
      )
      .select(
        'lineage_variant_primer_overlaps.ref_start', 'lineage_variant_primer_overlaps.ref_end',
        'lineage_variant_primer_overlaps.variant_type', 'lineage_variant_primer_overlaps.variant',
        'lineage_variant_primer_overlaps.ref', 'lineage_variant_primer_overlaps.frequency_pct',
        'lineage_variant_primer_overlaps.first_seen_date', 'lineage_variant_primer_overlaps.first_seen_lineage',
        'lineage_variant_primer_overlaps.first_seen_location', 'lineage_variant_primer_overlaps.last_seen_date',
        'lineage_variant_primer_overlaps.last_seen_lineage', 'lineage_variant_primer_overlaps.last_seen_location',
        'oligos.name AS oligo_name', 'oligos.sequence AS oligo_sequence',
        'oligos.strand AS oligo_strand', 'oligo_alignment_positions.ref_start AS oligo_start',
        'oligo_alignment_positions.ref_end AS oligo_end', 'primer_sets.name AS primer_set_name'
      )
      .order('lineage_variant_primer_overlaps.ref_start ASC, lineage_variant_primer_overlaps.ref_end ASC')
  end
end
