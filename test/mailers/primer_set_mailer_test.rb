# frozen_string_literal: true

require 'test_helper'

class PrimerSetMailerTest < ActionMailer::TestCase
  test 'updated_primer_set_email' do
    # Create the email and store it for further assertions
    email = PrimerSetMailer.updated_primer_set_email('me@example.com', primer_sets(:one))

    # Send the email, then test that it got queued
    assert_emails 1 do
      email.deliver_now
    end

    # Test the body of the sent email contains what we expect it to
    assert_equal ['me@example.com'], email.to
    assert_equal [ENV.fetch('ADMIN_EMAIL', nil)], email.from
    assert_equal 'Primer Set Created or Updated', email.subject
  end

  test 'primer_overlap_notification_email links each position to its variant summary' do
    user = users(:one)
    primer_set = primer_sets(:one)
    oligo = oligos(:one)
    notification = ProposedNotification.new(
      user:, primer_set:, oligo:, coordinate: 100, fraction_variant: 0.25,
      primer_set_subscription: PrimerSetSubscription.create!(user:, primer_set:),
      subscribed_geo_location: user.subscribed_geo_locations.create!(
        detailed_geo_location_alias: detailed_geo_location_aliases(:darwin)
      ),
      detailed_geo_location_alias: detailed_geo_location_aliases(:darwin)
    )
    notification.save!

    email = PrimerSetMailer.primer_overlap_notification_email(user.id, [notification])
    url = Rails.application.routes.url_helpers.organism_variant_position_url(
      organism_slug: primer_set.organism.slug, coordinate: 100, days: user.lookback_days, host: 'www.example.com'
    )

    assert_includes email.body.to_s, url
  end
end
