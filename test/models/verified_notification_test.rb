# frozen_string_literal: true

require 'test_helper'

class VerifiedNotificationTest < ActiveSupport::TestCase
  include ActionMailer::TestHelper

  test 'deliver! emails the user and records it as sent' do
    notification = VerifiedNotification.create!(user: users(:one), status: 'Unsent')

    assert_emails(1) { notification.deliver! }
    assert_equal 'Sent', notification.reload.status
  end

  test 'deliver! sends inactive users nothing and records the notification as skipped' do
    users(:one).update!(active: false)
    notification = VerifiedNotification.create!(user: users(:one), status: 'Unsent')

    assert_no_emails { notification.deliver! }
    assert_equal 'Skipped', notification.reload.status
  end

  test 'skip_unsent! marks pending notifications handled, except for the given addresses' do
    skipped = VerifiedNotification.create!(user: users(:one), status: 'Unsent')
    kept = VerifiedNotification.create!(user: users(:two), status: 'Unsent')

    VerifiedNotification.skip_unsent!(except_emails: [" #{users(:two).email.upcase}"])

    assert_equal 'Skipped', skipped.reload.status
    assert_equal 'Unsent', kept.reload.status
    assert_equal [kept], VerifiedNotification.find_or_create_verified_notifications!, 'only the kept one is sent'
  end

  test 'send_summary! emails each active administrator about handled notifications' do
    sent = VerifiedNotification.create!(user: users(:one), status: 'Sent')
    pending = VerifiedNotification.create!(user: users(:two), status: 'Unsent')

    assert_emails(1) { VerifiedNotification.send_summary!([sent, pending]) }

    users(:admin_user).update!(active: false)

    assert_no_emails { VerifiedNotification.send_summary!([sent]) }
  end

  test 'send_summary! sends nothing when no notification was handled' do
    pending = VerifiedNotification.create!(user: users(:two), status: 'Unsent')

    assert_no_emails { VerifiedNotification.send_summary!([]) }
    assert_no_emails { VerifiedNotification.send_summary!([pending]) }
  end
end
