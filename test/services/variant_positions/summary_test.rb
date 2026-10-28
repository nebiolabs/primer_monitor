# frozen_string_literal: true

require 'test_helper'

module VariantPositions
  class SummaryTest < ActiveSupport::TestCase
    setup do
      @organism = organisms(:sars_cov2)
      @taxon = organism_taxa(:sars_cov_2_taxon)
      @darwin = detailed_geo_locations(:darwin)
      @warumungu = detailed_geo_locations(:warumungu)
      @xbb = Lineage.create!(name: 'XBB.1', organism: @organism)
      @xbb_call = LineageCall.create!(taxon: 'XBB_ACC', lineage: @xbb)

      @snp_carrier = record('snp_carrier', @darwin, lineage_call: @xbb_call)
      @clean = record('clean', @darwin, lineage_call: @xbb_call)
      @deletion_carrier = record('deletion_carrier', @warumungu)
      @old_carrier = record('old_carrier', @darwin, date: 200.days.ago.to_date)

      site(@snp_carrier, ref_start: 100, ref_end: 101, variant_type: 'X', variant: 'A')
      site(@deletion_carrier, ref_start: 98, ref_end: 103, variant_type: 'D', variant: '5')
      site(@old_carrier, ref_start: 100, ref_end: 101, variant_type: 'X', variant: 'G')
    end

    def record(strain, location, lineage_call: nil, date: Date.current)
      FastaRecord.create!(strain:, detailed_geo_location: location, date_collected: date,
                          organism_taxon_id: @taxon.id, lineage_call:)
    end

    def site(fasta_record, **attrs)
      VariantSite.create!(fasta_record_id: fasta_record.id, organism_taxon_id: @taxon.id, ref: nil, **attrs)
    end

    def summarize(coordinate, days: 30)
      Summary.call(organism: @organism, coordinate:, days:)
    end

    test 'counts sequences carrying a variant at the position within the window' do
      summary = summarize(100)

      assert_equal 2, summary[:variant_sequences]
      assert_equal 3, summary[:sequences]
    end

    test 'breaks variant sequences down by location with the fraction of the location' do
      darwin, warumungu = summarize(100)[:locations]

      assert_equal 'Oceania/Australia/Northern Territory/Darwin', darwin.name
      assert_equal [1, 2], [darwin.variant_count, darwin.sequence_count]
      assert_in_delta 0.5, darwin.fraction
      assert_equal 'Oceania/Australia/Northern Territory/Warumungu', warumungu.name
      assert_in_delta 1.0, warumungu.fraction
    end

    test 'locations carry their coordinates for the map' do
      darwin = summarize(100)[:locations].first

      assert_in_delta detailed_geo_location_aliases(:darwin).latitude, darwin.latitude
      assert_in_delta detailed_geo_location_aliases(:darwin).longitude, darwin.longitude
    end

    test 'lists only named lineages that carry the variant' do
      lineage = summarize(100)[:lineages].first

      assert_equal 'XBB.1', lineage.name
      assert_equal [1, 2], [lineage.variant_count, lineage.sequence_count]
    end

    test 'a deletion covers every position in its span, the SNP ends before position 101' do
      assert_equal 1, summarize(101)[:variant_sequences]
      assert_equal 1, summarize(98)[:variant_sequences]
    end

    test 'no variants outside the position' do
      summary = summarize(500)

      assert_equal 0, summary[:variant_sequences]
      assert_empty summary[:locations]
      assert_empty summary[:lineages]
    end

    test 'an insertion counts only at its start position' do
      site(@clean, ref_start: 300, ref_end: 304, variant_type: 'I', variant: 'ACGT')

      assert_equal 1, summarize(300)[:variant_sequences]
      assert_equal 0, summarize(301)[:variant_sequences]
    end

    test 'a longer window includes older sequences' do
      assert_equal 3, summarize(100, days: 365)[:variant_sequences]
    end
  end
end
