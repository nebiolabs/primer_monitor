# frozen_string_literal: true

require 'test_helper'

class SubscriptionsTest < ActionDispatch::IntegrationTest
  setup do
    @user = users(:one)
    sign_in(@user)
    @subscription = PrimerSetSubscription.create!(user: @user, primer_set: primer_sets(:one), active: true)
    PrimerSetSubscription.create!(user: @user, primer_set: primer_sets(:two), active: false)
    darwin = detailed_geo_location_aliases(:darwin)
    location = SubscribedGeoLocation.create!(user: @user, detailed_geo_location_alias: darwin)
    ProposedNotification.create!(user: @user, primer_set: primer_sets(:one), oligo: oligos(:one), coordinate: 100,
                                 fraction_variant: 0.25, primer_set_subscription: @subscription,
                                 subscribed_geo_location: location, detailed_geo_location_alias: darwin,
                                 verified_notification: VerifiedNotification.create!(user: @user, status: 'Sent'))
  end

  test 'summaries cover active subscriptions only, with what was sent' do
    summary, *others = PrimerSetSubscription.summaries_for(@user)

    assert_empty others
    assert_equal primer_sets(:one), summary.primer_set
    assert_equal 1, summary.notifications_sent
    assert_not_nil summary.last_notified_at
    assert_nil summary.overlapping_variants, 'the overlap views are not populated in test'
  end

  test 'the account page lists active subscriptions and the linked sign-in' do
    @user.update!(provider: 'entra_id', uid: 'uid-1')

    get edit_user_registration_path

    assert_response :success
    assert_select '.subscriptions-table tbody tr', 1
    assert_select '.subscriptions-table a', text: primer_sets(:one).name
    assert_match 'Microsoft', response.body
  end

  test 'the primer set page shows the notification history' do
    get primer_set_path(primer_sets(:one))

    assert_select '.notification-history tbody tr', 1
    assert_select '.notification-history td', text: 'Sent'
  end

  test 'subscribing to a primer set that does not exist is a 404 and changes nothing' do
    post primer_set_subscriptions_path(primer_set_id: 0)

    assert_response :not_found
    assert_not @user.reload.send_primer_updates?
  end

  test 'unsubscribing from the account page returns there' do
    delete primer_set_subscription_path(@subscription), headers: { 'HTTP_REFERER' => edit_user_registration_url }

    assert_redirected_to edit_user_registration_url
    assert_not @subscription.reload.active?
  end

  test 'the account page and password reset require signing in' do
    sign_out(@user)

    get edit_user_registration_path

    assert_redirected_to new_user_session_path
    assert_no_emails { post user_password_reset_path }
  end

  test 'a signed-in user can ask for a password reset email' do
    assert_emails(1) { post user_password_reset_path }

    assert_redirected_to edit_user_registration_path
  end

  test 'location options list each location once' do
    world = detailed_geo_location_aliases(:darwin)
    DetailedGeoLocationAlias.stubs(:subscribable).returns([world, world])

    assert_equal [[world.name, world.id]], DetailedGeoLocationAlias.subscribable_options
  end
end
