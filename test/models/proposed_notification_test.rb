# frozen_string_literal: true

require 'test_helper'

class ProposedNotificationTest < ActiveSupport::TestCase
  setup do
    location = detailed_geo_locations(:darwin)
    fasta_ids = Array.new(20) do |i|
      FastaRecord.create!(strain: "test/#{i}", detailed_geo_location: location, date_collected: Date.current,
                          organism_taxon_id: organism_taxa(:sars_cov_2_taxon).id).id
    end
    fasta_ids.first(5).each do |id|
      VariantSite.create!(fasta_record_id: id, ref_start: 100, ref_end: 101, variant_type: 'X', variant: 'A')
    end
    @location_sub = users(:one).subscribed_geo_locations
                               .create!(detailed_geo_location_alias: detailed_geo_location_aliases(:darwin))
    @primer_set_sub = PrimerSetSubscription.create!(user: users(:one), primer_set: primer_sets(:one))

    %w[variant_overlaps oligo_variant_overlaps identify_primers_for_notifications].each do |view|
      ActiveRecord::Base.connection.execute("REFRESH MATERIALIZED VIEW #{view}")
    end
  end

  def user_one_notifications
    ProposedNotification.new_proposed_notifications.select { |pn| pn.user_id == users(:one).id }
  end

  test 'builds a notification for a current subscription' do
    assert_not_empty user_one_notifications
  end

  # identify_primers_for_notifications is refreshed nightly; subscriptions can change before the morning send
  test 'skips a location dropped since the view refresh' do
    @location_sub.destroy!

    assert_empty user_one_notifications
  end

  test 'skips a primer set dropped since the view refresh' do
    @primer_set_sub.destroy!

    assert_empty user_one_notifications
  end
end
