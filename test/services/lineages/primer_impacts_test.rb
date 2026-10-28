# frozen_string_literal: true

require 'test_helper'

module Lineages
  class PrimerImpactsTest < ActiveSupport::TestCase
    setup do
      @organism = organisms(:sars_cov2)
      @taxon = organism_taxa(:sars_cov_2_taxon)
      primer_sets(:one).update!(status: 'complete')
      @lineage = Lineage.create!(name: 'XFG.TEST', organism: @organism)
      @call = LineageCall.create!(taxon: 'XFG_ACC', lineage: @lineage)
    end

    def sequences(count, variant_at: nil)
      count.times do |i|
        record = FastaRecord.create!(strain: "xfg/#{variant_at}/#{i}", organism_taxon_id: @taxon.id,
                                     detailed_geo_location: detailed_geo_locations(:darwin),
                                     date_collected: Date.current, lineage_call: @call)
        next if variant_at.nil?

        VariantSite.create!(fasta_record_id: record.id, organism_taxon_id: @taxon.id, ref: nil,
                            ref_start: 110, ref_end: 111, variant_type: 'X', variant: 'T')
      end
    end

    def impacts
      PrimerImpacts.new(organism: @organism, lineage: @lineage, since: Date.current - 30)
    end

    test 'a variant the lineage carries often enough scores the primers it overlaps' do
      sequences(12, variant_at: 'snp')
      ActiveRecord::Base.connection.execute('REFRESH MATERIALIZED VIEW lineage_info')

      primer_set = impacts.set(primer_sets(:one).id)

      # the SNP at 110 scores 3 on the primer, which has no recorded strand and reads as reverse: log2(2 ** 3) = 3
      assert_in_delta 3.0, primer_set.impact, 0.0001
      assert_in_delta 3.0, primer_set.primers.find { |p| p.ref_start == 95 }.impact, 0.0001
    end

    test 'a variant carried by fewer than 10 sequences does not count' do
      sequences(12, variant_at: 'rare')
      FastaRecord.where("strain LIKE 'xfg/rare/%'").order(:id).limit(4).each do |record|
        VariantSite.where(fasta_record_id: record.id).delete_all
      end

      assert_nil impacts.set(primer_sets(:one).id).impact
    end

    test 'the BED track names each primer with its impact and gives the set strand and score' do
      sequences(12, variant_at: 'bed')
      bed = impacts.bed(impacts.set(primer_sets(:one).id))

      fields = bed.split("\t")

      assert_equal @taxon.reference_accession, fields[0]
      assert_equal [95, 115], [fields[1].to_i, fields[2].to_i]
      assert_equal 3, fields[4].to_i
    end
  end
end
