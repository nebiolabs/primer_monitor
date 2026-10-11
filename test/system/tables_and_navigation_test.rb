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

  test 'primer sets can be filtered by amplification method' do
    primer_sets(:one).update!(status: :complete) # logged-out visitors only see complete sets
    primer_sets(:two).update!(status: :complete, amplification_method: amplification_methods(:tiled_amplicon))
    visit organism_primer_sets_path(organisms(:sars_cov2))

    assert_equal ['Primer Set', 'Amp. Method'], all('thead th .dt-column-title').first(2).map(&:text)

    find('thead th', text: 'Amp. Method').find('.dtcc-button_dropdown').click
    within('.dtcc-dropdown') { click_button 'Tiled amplicon' }

    assert_selector 'tbody tr', count: 1
    assert_selector 'tbody tr', text: 'CDC'
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
