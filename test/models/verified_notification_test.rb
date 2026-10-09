# frozen_string_literal: true

require 'test_helper'

class VerifiedNotificationTest < ActiveSupport::TestCase
  test 'skip_unsent! marks pending notifications handled, except for the given addresses' do
    skipped = VerifiedNotification.create!(user: users(:one), status: 'Unsent')
    kept = VerifiedNotification.create!(user: users(:two), status: 'Unsent')

    VerifiedNotification.skip_unsent!(except_emails: [" #{users(:two).email.upcase}"])

    assert_equal 'Skipped', skipped.reload.status
    assert_equal 'Unsent', kept.reload.status
    assert_equal [kept], VerifiedNotification.find_or_create_verified_notifications!, 'only the kept one is sent'
  end
end
