# frozen_string_literal: true

module Lineages
  # How a lineage's variants impact the organism's complete primer sets. A variant counts when the lineage carries it
  # at the classifier's threshold: at least 10 sequences and at least 1% of the lineage's sequences in the window. Each
  # counted variant is scored against the primers it overlaps (see VariantPositions::Disruption). A primer's impact is
  # the log2 of the sum of those contributions; a primer set's impact is that of its most affected primer.
  class PrimerImpacts
    MIN_SEQUENCES = 10
    MIN_FRACTION = 0.01
    CACHE_TTL = 1.hour

    Primer = Struct.new(:oligo_id, :set_id, :set_name, :name, :ref_start, :ref_end, :strand, :impact,
                        keyword_init: true)
    PrimerSet = Struct.new(:id, :name, :impact, :primers, keyword_init: true)

    def initialize(organism:, lineage:, since:)
      @organism = organism
      @lineage = lineage
      @since = since
      @taxon = organism.organism_taxa.first
    end

    # every complete primer set of the organism, most impacted first (sets without an impact last)
    def sets
      @sets ||= Rails.cache.fetch(['lineage_primer_impacts', lineage.id, taxon&.id, since], expires_in: CACHE_TTL) do
        build_sets
      end
    end

    def set(id)
      sets.find { |primer_set| primer_set.id == id }
    end

    # BED lines of one primer set's primers across the reference, named with their impact
    def bed(primer_set)
      primer_set.primers.map do |primer|
        impact = primer.impact ? format('%.1f', primer.impact) : 'none'
        [taxon.reference_accession, primer.ref_start, primer.ref_end, "#{primer.name}, impact #{impact}",
         (primer.impact || 0).round, primer.strand].join("\t").tr(' ', '_')
      end.join("\n")
    end

    private

    attr_reader :organism, :lineage, :since, :taxon

    def build_sets
      return [] if taxon.nil?

      primers = complete_primers(VariantSite.load_reference_sequence(taxon))
      impacts = primer_impacts(primers)
      sets = primers.group_by(&:set_id).map do |set_id, set_primers|
        set_primers.each { |primer| primer.impact = impacts[primer.oligo_id] }
        PrimerSet.new(id: set_id, name: set_primers.first.set_name, impact: set_primers.filter_map(&:impact).max,
                      primers: set_primers.sort_by(&:ref_start))
      end
      sets.sort_by { |primer_set| -(primer_set.impact || -Float::INFINITY) }
    end

    # the complete primer sets' primers, with strands settled against the reference
    def complete_primers(reference)
      positions = OligoAlignmentPosition.where(organism_taxon: taxon)
                                        .joins(oligo: :primer_set)
                                        .where(primer_sets: { organism_id: organism.id, status: 'complete' })
                                        .includes(oligo: :primer_set)
      positions.map do |position|
        oligo = position.oligo
        Primer.new(oligo_id: oligo.id, set_id: oligo.primer_set_id, set_name: oligo.primer_set.name,
                   name: oligo.label, ref_start: position.ref_start, ref_end: position.ref_end,
                   strand: VariantPositions::Strand.of(oligo, position.ref_start, position.ref_end, reference))
      end
    end

    # oligo id => impact, for the primers the lineage's frequent variants affect
    def primer_impacts(primers)
      scores = Hash.new { |hash, key| hash[key] = [] }
      frequent_variants.each do |start, stop, type|
        span = start...[stop, start + 1].max
        primers.each do |primer|
          next unless primer.ref_start < span.end && primer.ref_end > span.begin

          score = VariantPositions::Disruption.primer_score(variant_type: type, variant_span: span,
                                                            primer_span: primer.ref_start...primer.ref_end,
                                                            strand: primer.strand)
          scores[primer.oligo_id] << score if score
        end
      end
      scores.transform_values { |exponents| VariantPositions::Disruption.log2_sum(exponents) }
    end

    # [ref_start, ref_end, variant_type] of each distinct variant the lineage carries at the threshold
    def frequent_variants
      sql = ActiveRecord::Base.sanitize_sql_array([<<~SQL.squish, lineage.id, taxon.id, since, threshold])
        SELECT vs.ref_start, vs.ref_end, vs.variant_type
        FROM variant_sites vs
        JOIN fasta_records fr ON fr.id = vs.fasta_record_id
        JOIN lineage_calls lc ON lc.id = fr.lineage_call_id
        WHERE lc.lineage_id = ? AND vs.organism_taxon_id = ? AND fr.date_collected >= ?
          AND vs.variant_type IN ('X', 'D', 'I') AND vs.variant NOT LIKE '%N%'
        GROUP BY vs.ref_start, vs.ref_end, vs.variant_type, vs.variant
        HAVING count(*) >= ?
      SQL
      ActiveRecord::Base.connection.select_rows(sql).map { |start, stop, type| [start.to_i, stop.to_i, type] }
    end

    def threshold
      [MIN_SEQUENCES, (lineage_sequence_count * MIN_FRACTION).ceil].max
    end

    def lineage_sequence_count
      @lineage_sequence_count ||= FastaRecord.joins(lineage_call: :lineage)
                                             .where(lineage_calls: { lineage_id: lineage.id },
                                                    organism_taxon_id: taxon.id, date_collected: since..)
                                             .count
    end
  end
end
