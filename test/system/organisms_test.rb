# frozen_string_literal: true

require 'application_system_test_case'

class OrganismsTest < ApplicationSystemTestCase
  setup do
    @organism = organisms(:one)
  end

  test 'visiting the index' do
    visit organisms_url
    assert_selector 'h1', text: 'Organisms'
  end

  test 'creating a Organism' do
    visit organisms_url
    click_on 'New Organism'

    fill_in 'Name', with: @organism.name
    fill_in 'Ncbi taxon', with: @organism.ncbi_taxon_id
    click_on 'Create Organism'

    assert_text 'Organism was successfully created'
    click_on 'Back'
  end

  test 'updating a Organism' do
    visit organisms_url
    click_on 'Edit', match: :first

    fill_in 'Name', with: @organism.name
    fill_in 'Ncbi taxon', with: @organism.ncbi_taxon_id
    click_on 'Update Organism'

    assert_text 'Organism was successfully updated'
    click_on 'Back'
  end

  test 'destroying a Organism' do
    visit organisms_url
    page.accept_confirm do
      click_on 'Destroy', match: :first
    end

    assert_text 'Organism was successfully destroyed'
  end
end
