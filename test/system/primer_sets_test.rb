require "application_system_test_case"

class PrimerSetsTest < ApplicationSystemTestCase
  setup do
    @primer_set = primer_sets(:one)
  end

  test "visiting the index" do
    visit primer_sets_url
    assert_selector "h1", text: "Primer Sets"
  end

  test "creating a Primer set" do
    visit primer_sets_url
    click_on "New Primer Set"

    fill_in "Creator", with: @primer_set.creator_id
    fill_in "Name", with: @primer_set.name
    click_on "Create Primer set"

    assert_text "Primer set was successfully created"
    click_on "Back"
  end

  test "updating a Primer set" do
    visit primer_sets_url
    click_on "Edit", match: :first

    fill_in "Creator", with: @primer_set.creator_id
    fill_in "Name", with: @primer_set.name
    click_on "Update Primer set"

    assert_text "Primer set was successfully updated"
    click_on "Back"
  end

  test "destroying a Primer set" do
    visit primer_sets_url
    page.accept_confirm do
      click_on "Destroy", match: :first
    end

    assert_text "Primer set was successfully destroyed"
  end
end
