# frozen_string_literal: true

require 'test_helper'

module VariantPositions
  class TrackTest < ActiveSupport::TestCase
    setup do
      @organism = organisms(:sars_cov2)
      @taxon = organism_taxa(:sars_cov_2_taxon)
      @track = Track.new(organism: @organism, coordinate: 100)
    end

    test 'primer track lists the primers covering the position' do
      lines = @track.primer_bed.split("\n")

      assert_equal 1, lines.size

      chrom, start, stop, name = lines.first.split("\t")

      assert_equal @taxon.reference_accession, chrom
      assert_equal [95, 115], [start.to_i, stop.to_i]
      assert_includes name, oligos(:one).label
    end

    test 'variant track marks the variants at the position with how many sequences carry them' do
      2.times do |i|
        record = FastaRecord.create!(strain: "track/#{i}", detailed_geo_location: detailed_geo_locations(:darwin),
                                     date_collected: Date.current, organism_taxon_id: @taxon.id)
        VariantSite.create!(fasta_record_id: record.id, organism_taxon_id: @taxon.id, ref: nil,
                            ref_start: 100, ref_end: 101, variant_type: 'X', variant: 'A')
      end
      other = FastaRecord.create!(strain: 'track/other', detailed_geo_location: detailed_geo_locations(:darwin),
                                  date_collected: Date.current, organism_taxon_id: @taxon.id)
      VariantSite.create!(fasta_record_id: other.id, organism_taxon_id: @taxon.id, ref: nil,
                          ref_start: 300, ref_end: 301, variant_type: 'X', variant: 'G')

      assert_equal "#{@taxon.reference_accession}\t100\t101\tX:A (2 seqs)\t0\t+", @track.variant_bed
    end

    test 'locus frames the covering primers with padding, in 1-based coordinates' do
      assert_equal "#{@taxon.reference_accession}:71-140", @track.locus
    end

    test 'locus falls back to a window around the position when no primer covers it' do
      track = Track.new(organism: @organism, coordinate: 5000)

      assert_equal "#{@taxon.reference_accession}:4976-5026", track.locus
    end
  end
end
