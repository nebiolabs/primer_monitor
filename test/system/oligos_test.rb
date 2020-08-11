require "application_system_test_case"

class OligosTest < ApplicationSystemTestCase
  setup do
    @oligo = oligos(:one)
  end

  test "visiting the index" do
    visit oligos_url
    assert_selector "h1", text: "Oligos"
  end

  test "creating a Oligo" do
    visit oligos_url
    click_on "New Oligo"

    fill_in "Name", with: @oligo.name
    fill_in "Sequence", with: @oligo.sequence
    click_on "Create Oligo"

    assert_text "Oligo was successfully created"
    click_on "Back"
  end

  test "updating a Oligo" do
    visit oligos_url
    click_on "Edit", match: :first

    fill_in "Name", with: @oligo.name
    fill_in "Sequence", with: @oligo.sequence
    click_on "Update Oligo"

    assert_text "Oligo was successfully updated"
    click_on "Back"
  end

  test "destroying a Oligo" do
    visit oligos_url
    page.accept_confirm do
      click_on "Destroy", match: :first
    end

    assert_text "Oligo was successfully destroyed"
  end
end
