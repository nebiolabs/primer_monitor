# frozen_string_literal: true

require 'test_helper'

# The notification fraction: sequences with a variant at a primer position over all recent sequences there.
class IdentifyPrimersForNotificationTest < ActiveSupport::TestCase
  test 'overlapping location subscriptions and other subscribers count each sequence once' do
    location = detailed_geo_locations(:darwin)
    fasta_ids = Array.new(20) do |i|
      FastaRecord.create!(strain: "test/#{i}", detailed_geo_location: location, date_collected: Date.current,
                          organism_taxon_id: organism_taxa(:sars_cov_2_taxon).id).id
    end
    # 5 of the 20 have a SNP inside oligo one (aligned at 95-115)
    fasta_ids.first(5).each do |id|
      VariantSite.create!(fasta_record_id: id, ref_start: 100, ref_end: 101, variant_type: 'X', variant: 'A')
    end
    # World and Darwin both cover the Darwin location
    users(:one).subscribed_geo_locations.create!(detailed_geo_location_alias: detailed_geo_location_aliases(:world))
    users(:one).subscribed_geo_locations.create!(detailed_geo_location_alias: detailed_geo_location_aliases(:darwin))
    users(:two).subscribed_geo_locations.create!(detailed_geo_location_alias: detailed_geo_location_aliases(:darwin))
    PrimerSetSubscription.create!(user: users(:one), primer_set: primer_sets(:one))

    %w[variant_overlaps oligo_variant_overlaps identify_primers_for_notifications].each do |view|
      ActiveRecord::Base.connection.execute("REFRESH MATERIALIZED VIEW #{view}")
    end

    row = IdentifyPrimersForNotification.where(user_id: users(:one).id, oligo_id: oligos(:one).id).sole

    assert_equal [100, 5, 20], [row.coords.to_i, row.variant_count, row.records_count]
    assert_in_delta 0.25, row.fraction_variant
  end

  test 'variant_overlaps keeps only the 168 days before the newest collection date' do
    taxon_id = organism_taxa(:sars_cov_2_taxon).id
    newest = Date.current # the window is counted back from this newest sequence, not from today
    dates = { newest:, recent: newest - 168 + 1, old: newest - 168 }
    ids = dates.transform_values do |date|
      FastaRecord.create!(strain: "test/#{date}", detailed_geo_location: detailed_geo_locations(:darwin),
                          date_collected: date, organism_taxon_id: taxon_id).id
    end
    ids.each_value do |id|
      VariantSite.create!(fasta_record_id: id, ref_start: 100, ref_end: 101, variant_type: 'X', variant: 'A')
    end

    ActiveRecord::Base.connection.execute('REFRESH MATERIALIZED VIEW variant_overlaps')
    variant_ids = VariantSite.where(fasta_record_id: ids.values).pluck(:fasta_record_id, :id).to_h
    in_view = ActiveRecord::Base.connection.select_values('SELECT DISTINCT variant_id FROM variant_overlaps')

    assert_includes in_view, variant_ids[ids[:recent]]
    assert_not_includes in_view, variant_ids[ids[:old]]
  end
end
