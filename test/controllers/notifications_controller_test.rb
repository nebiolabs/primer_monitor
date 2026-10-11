# frozen_string_literal: true

require 'test_helper'

class NotificationsControllerTest < ActionDispatch::IntegrationTest
  setup do
    @sent = VerifiedNotification.create!(user: users(:one), status: 'Sent')
    @skipped = VerifiedNotification.create!(user: users(:two), status: 'Skipped')
    sign_in(users(:admin_user))
  end

  test 'admins see who was emailed' do
    get notifications_url

    assert_response :success
    assert_select 'td', text: users(:one).email
    assert_select 'td', text: users(:two).email
  end

  test 'filters by status and user' do
    get notifications_url(status: 'Skipped')

    assert_select 'td', text: users(:two).email
    assert_select 'td', text: users(:one).email, count: 0

    get notifications_url(q: 'frank')

    assert_select 'td', text: users(:one).email
    assert_select 'td', text: users(:two).email, count: 0
  end

  test 'filters by date and ignores malformed dates' do
    get notifications_url(from: 1.week.from_now.to_date, to: 2.weeks.from_now.to_date)

    assert_select 'td', text: users(:one).email, count: 0

    get notifications_url(from: 'not-a-date')

    assert_select 'td', text: users(:one).email
  end

  test 'admins see an email in detail' do
    get notification_url(@sent)

    assert_response :success
    assert_select 'p', text: /#{users(:one).email}/
  end

  test 'other users cannot see the dashboard' do
    sign_out(users(:admin_user))
    sign_in(users(:one))
    get notifications_url

    assert_redirected_to root_url

    get notification_url(@sent)

    assert_redirected_to root_url
  end

  test 'visitors who are not logged in are sent to log in' do
    sign_out(users(:admin_user))
    get notifications_url

    assert_redirected_to new_user_session_url

    get notification_url(@sent)

    assert_redirected_to new_user_session_url
  end
end
