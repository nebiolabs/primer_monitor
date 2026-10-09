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

    %w[variant_overlaps counts time_counts oligo_variant_overlaps identify_primers_for_notifications].each do |view|
      ActiveRecord::Base.connection.execute("REFRESH MATERIALIZED VIEW #{view}")
    end

    row = IdentifyPrimersForNotification.where(user_id: users(:one).id, oligo_id: oligos(:one).id).sole

    assert_equal [100, 5, 20], [row.coords.to_i, row.variant_count, row.records_count]
    assert_in_delta 0.25, row.fraction_variant
  end
end
