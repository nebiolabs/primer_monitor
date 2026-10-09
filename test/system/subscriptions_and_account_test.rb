# frozen_string_literal: true

require 'application_system_test_case'

# Subscribing (Turbo forms, no rails-ujs), the account page's location picker (Tom Select) and signing out.
class SubscriptionsAndAccountTest < ApplicationSystemTestCase
  setup { sign_in users(:one) }

  test 'subscribing and unsubscribing on a primer set page' do
    visit primer_set_path(primer_sets(:two))
    click_on 'Subscribe'

    assert_button 'Unsubscribe'

    click_on 'Unsubscribe'

    assert_button 'Subscribe'
  end

  test 'the location picker lists each location once' do
    visit edit_user_registration_path

    options = page.evaluate_script(<<~JS)
      Object.values(document.getElementById('user_subscribed_detailed_geo_location_alias_ids').tomselect.options)
        .map(option => option.text)
    JS

    assert_not_empty options
    assert_equal options.uniq, options
  end

  test 'saving the user form updates the user' do
    visit edit_user_path(users(:one))
    fill_in 'user_first', with: 'Francis'
    click_on 'Save'

    assert_text 'User was successfully updated.'
    assert_equal 'Francis', users(:one).reload.first
  end

  test 'logging out' do
    visit root_path
    click_on 'Log out'

    assert_link 'Log in'
  end
end
