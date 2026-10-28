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
      lines = @track.primer_bed.split("\n")

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

      score = @track.primer_bed.split("\t")[4]

      assert_equal '1', score
    end

    test 'primer track gives each primer its strand' do
      oligos(:one).update!(strand: '-')

      assert_equal '-', @track.primer_bed.split("\t").last
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

      assert_equal "#{@taxon.reference_accession}#{tab}100#{tab}101#{tab}X:A_(50.0%)#{tab}2#{tab}.",
                   @track.variant_bed(0)
      assert_equal "#{@taxon.reference_accession}#{tab}100#{tab}101#{tab}X:G_(25.0%)#{tab}1#{tab}.",
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
