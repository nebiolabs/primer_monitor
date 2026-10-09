# frozen_string_literal: true

require 'application_system_test_case'

# Creating an account with a password: sign up, confirm by email, then log in.
class SignUpTest < ApplicationSystemTestCase
  test 'signing up, confirming the email address, and logging in' do
    visit new_user_registration_path
    fill_in 'user_first', with: 'Rosalind'
    fill_in 'user_last', with: 'Franklin'
    fill_in 'user_email', with: 'rosalind@example.org'
    fill_in 'user_password', with: 'photo51-password'
    fill_in 'user_password_confirmation', with: 'photo51-password'

    assert_difference('User.count') do
      click_on 'Sign up'

      assert_text 'confirmation link'
    end

    log_in 'rosalind@example.org', 'photo51-password'

    assert_text 'You have to confirm your email address'

    confirmation = ActionMailer::Base.deliveries.last

    assert_equal ['rosalind@example.org'], confirmation.to
    visit URI(confirmation.body.encoded[%r{https?://\S+confirmation_token=[\w-]+}]).request_uri

    assert_text 'successfully confirmed'

    log_in 'rosalind@example.org', 'photo51-password'

    assert_button 'Log out'
  end

  def log_in(email, password)
    visit new_user_session_path
    fill_in 'user_email', with: email
    fill_in 'user_password', with: password
    click_button 'Log in'
  end
end
