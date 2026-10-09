# frozen_string_literal: true

require 'application_system_test_case'

# DataTables (search, Copy/Excel export) and the mobile navigation menu.
class TablesAndNavigationTest < ApplicationSystemTestCase
  test 'tables get search and Copy/Excel export' do
    visit organism_primer_sets_path(organisms(:sars_cov2))

    assert_selector '.dt-search input'
    assert_button 'Copy'
    assert_button 'Excel'
  end

  test 'the burger opens and closes the mobile menu' do
    page.driver.browser.manage.window.resize_to(600, 900)
    visit root_path
    find('.navbar-burger').click

    assert_selector '.navbar-menu.is-active'

    find('.navbar-burger').click

    assert_no_selector '.navbar-menu.is-active'
  ensure
    page.driver.browser.manage.window.resize_to(1400, 1400)
  end
end
