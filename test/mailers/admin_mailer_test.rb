# frozen_string_literal: true

require 'test_helper'

class AdminMailerTest < ActionMailer::TestCase
  test 'notification_summary_email lists who was emailed and links to the dashboard for that day' do
    sent = VerifiedNotification.create!(user: users(:one), status: 'Sent')
    skipped = VerifiedNotification.create!(user: users(:two), status: 'Skipped')
    date = Date.new(2026, 10, 11)

    email = AdminMailer.notification_summary_email(users(:admin_user), [sent, skipped], date)

    assert_emails(1) { email.deliver_now }
    assert_equal [users(:admin_user).email], email.to
    assert_equal 'Primer Monitor: 1 notification email sent 2026-10-11', email.subject
    assert_includes email.body.to_s, users(:one).email
    assert_includes email.body.to_s, users(:two).email
    assert_includes email.body.to_s, 'notifications?from=2026-10-11&amp;to=2026-10-11'
  end
end
