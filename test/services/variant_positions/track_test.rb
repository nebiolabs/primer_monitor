# frozen_string_literal: true

require 'test_helper'

module VariantPositions
  class TrackTest < ActiveSupport::TestCase
    setup do
      @organism = organisms(:sars_cov2)
      @taxon = organism_taxa(:sars_cov_2_taxon)
      @track = Track.new(organism: @organism, coordinate: 100, days: 30)
      refresh_overlap_views
    end

    def refresh_overlap_views
      %w[variant_overlaps oligo_variant_overlaps].each do |view|
        ActiveRecord::Base.connection.execute("REFRESH MATERIALIZED VIEW #{view}")
      end
    end

    test 'primer track lists the primers overlapping the position' do
      lines = @track.primer_bed(0).split("\n")

      assert_equal 1, lines.size

      chrom, start, stop, name, _score, strand = lines.first.split("\t")

      assert_equal @taxon.reference_accession, chrom
      assert_equal [95, 115], [start.to_i, stop.to_i]
      assert_includes name, oligos(:one).label.tr(' ', '_')
      assert_equal '.', strand
    end

    test 'primer score counts the sequences in the window with a variant inside the primer' do
      record = FastaRecord.create!(strain: 'track/primer', organism_taxon_id: @taxon.id,
                                   detailed_geo_location: detailed_geo_locations(:darwin), date_collected: Date.current)
      VariantSite.create!(fasta_record_id: record.id, organism_taxon_id: @taxon.id, ref: nil,
                          ref_start: 110, ref_end: 111, variant_type: 'X', variant: 'T')
      VariantSite.create!(fasta_record_id: record.id, organism_taxon_id: @taxon.id, ref: nil,
                          ref_start: 300, ref_end: 301, variant_type: 'X', variant: 'T')
      refresh_overlap_views

      score = @track.primer_bed(0).split("\t")[4]

      assert_equal '3', score # the SNP scores 3 at the primer's 3' end, so its impact is log2(2 ** 3) = 3
    end

    test 'an allele is coloured by its disruption to the overlapping primers' do
      record = FastaRecord.create!(strain: 'track/disruption', organism_taxon_id: @taxon.id,
                                   detailed_geo_location: detailed_geo_locations(:darwin), date_collected: Date.current)
      VariantSite.create!(fasta_record_id: record.id, organism_taxon_id: @taxon.id, ref: nil,
                          ref_start: 100, ref_end: 101, variant_type: 'X', variant: 'A')
      allele = @track.alleles.first

      assert_equal 4, @track.disruption(allele) # the primer's strand is unset, so it is read as reverse
      assert_equal '#f3b758', @track.allele_color(allele) # a fifth of the way from yellow to red
      assert_match(/disruption 4/, @track.allele_name(allele))
    end

    test 'primer impact combines its variants as log2 of the sum of 2 ** disruption, and colours the primer' do
      assert_nil @track.primers.first.impact
      assert_equal Disruption::NO_IMPACT_COLOR, @track.primer_color(@track.primers.first)

      [110, 100].each do |start|
        record = FastaRecord.create!(strain: "track/impact#{start}", organism_taxon_id: @taxon.id,
                                     detailed_geo_location: detailed_geo_locations(:darwin),
                                     date_collected: Date.current)
        VariantSite.create!(fasta_record_id: record.id, organism_taxon_id: @taxon.id, ref: nil,
                            ref_start: start, ref_end: start + 1, variant_type: 'X', variant: 'T')
      end
      refresh_overlap_views

      track = Track.new(organism: @organism, coordinate: 100, days: 30) # a fresh track, after the variants exist
      primer = track.primers.first
      # on the primer's unset (reverse) strand a SNP at 110 scores 3 and one at 100 scores 4
      assert_in_delta Math.log2((2**3) + (2**4)), primer.impact, 0.0001
      assert_equal Disruption.color(primer.impact), track.primer_color(primer)
    end

    test 'primer track gives each primer its strand' do
      oligos(:one).update!(strand: '-')

      assert_equal '-', @track.primer_bed(0).split("\t").last
    end

    test 'variant tracks give each allele at the position its share of the window, most frequent first' do
      { 'A' => 2, 'G' => 1 }.each do |variant, count|
        count.times do |i|
          record = FastaRecord.create!(strain: "track/#{variant}#{i}", organism_taxon_id: @taxon.id,
                                       detailed_geo_location: detailed_geo_locations(:darwin),
                                       date_collected: Date.current)
          VariantSite.create!(fasta_record_id: record.id, organism_taxon_id: @taxon.id, ref: nil,
                              ref_start: 100, ref_end: 101, variant_type: 'X', variant:)
        end
      end
      FastaRecord.create!(strain: 'track/clean', detailed_geo_location: detailed_geo_locations(:darwin),
                          date_collected: Date.current, organism_taxon_id: @taxon.id)

      tab = "\t"

      assert_equal "#{@taxon.reference_accession}#{tab}100#{tab}101#{tab}X:A_(50.0%),_disruption_4#{tab}2#{tab}.",
                   @track.variant_bed(0)
      assert_equal "#{@taxon.reference_accession}#{tab}100#{tab}101#{tab}X:G_(25.0%),_disruption_4#{tab}1#{tab}.",
                   @track.variant_bed(1)
      assert_nil @track.variant_bed(2)
    end

    test 'locus frames the overlapping primers with padding, in 1-based coordinates' do
      assert_equal "#{@taxon.reference_accession}:71-140", @track.locus
    end

    test 'locus falls back to a window around the position when no primer overlaps it' do
      track = Track.new(organism: @organism, coordinate: 5000, days: 30)

      assert_equal "#{@taxon.reference_accession}:4976-5026", track.locus
    end
  end
end
