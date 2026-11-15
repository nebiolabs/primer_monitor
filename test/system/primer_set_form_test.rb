# frozen_string_literal: true

require 'application_system_test_case'

# The primer set form's oligo table: adding and removing rows (nested_fields.js), and sorting, bulk editing, pasting
# and FASTA loading (oligo_editor.js).
class PrimerSetFormTest < ApplicationSystemTestCase
  setup { sign_in users(:one) }

  def oligo_rows = all('#samples .nested-fields')

  test 'oligo rows can be added and removed' do
    visit new_primer_set_path
    click_on 'Add Oligo'
    click_on 'Add Oligo'

    assert_equal 2, oligo_rows.size

    oligo_rows.last.find('[data-nested-remove]').click

    assert_equal 1, oligo_rows.size
  end

  test 'uploading a FASTA file adds one oligo row per sequence, and the set saves with them' do
    fasta = Rails.root.join('tmp/system_test_primers.fasta')
    File.write(fasta, ">Fwd_1 forward primer\nACGTACGTACGTACGTACGT\n>Rev_1 reverse primer\nTTGGCCAATTGGCCAATTGG\n")

    visit new_primer_set_path
    fill_in 'primer_set_name', with: 'Uploaded set'
    select amplification_methods(:qPCR).name, from: 'primer_set_amplification_method_id'
    attach_file 'fasta_upload', fasta, make_visible: true

    assert_selector '#samples .nested-fields', count: 2
    assert_equal %w[ACGTACGTACGTACGTACGT TTGGCCAATTGGCCAATTGG], all('#samples input[id$="_sequence"]').map(&:value)

    assert_difference('Oligo.count', 2) do
      click_on 'Create Primer set'

      assert_text 'successfully added'
    end
  ensure
    FileUtils.rm_f(fasta)
  end

  def oligo_names = all('#samples input[data-col="name"]').map(&:value)

  def edit_set_with(names)
    primer_set = primer_sets(:one)
    primer_set.oligos.destroy_all
    names.each { |name| primer_set.oligos.create!(name:, sequence: 'ACGTACGTAC') }
    visit edit_primer_set_path(primer_set)
    primer_set
  end

  test 'oligos sort in natural order and the header toggles the direction' do
    edit_set_with(%w[amp_10_LEFT amp_2_LEFT amp_1_LEFT])

    assert_equal %w[amp_1_LEFT amp_2_LEFT amp_10_LEFT], oligo_names

    find('th', text: 'Name').click

    assert_equal %w[amp_10_LEFT amp_2_LEFT amp_1_LEFT], oligo_names
  end

  test 'a bulk category change applies to the checked rows that the filter shows' do
    primer_set = edit_set_with(%w[a_1_LEFT a_1_RIGHT b_1_LEFT])
    find('[data-oligo-filter]').fill_in(with: 'a_1')
    find('[data-select-all]').check
    find('[data-bulk-category]').select('Probe')
    click_on 'Set category'
    click_on 'Update Primer set'

    assert_text 'successfully updated'
    assert_equal({ 'a_1_LEFT' => 'Probe', 'a_1_RIGHT' => 'Probe', 'b_1_LEFT' => 'Forward' },
                 primer_set.oligos.reload.to_h { |o| [o.name, o.category] })
  end

  test 'pasting several lines into a cell fills the rows below and adds rows as needed' do
    edit_set_with(%w[x_1_LEFT])
    first('#samples input[data-col="sequence"]').click
    page.execute_script(<<~'JS')
      const input = document.querySelector('#samples input[data-col="sequence"]');
      const data = new DataTransfer();
      data.setData('text/plain', 'AAAACCCC\nGGGGTTTT\n');
      input.dispatchEvent(new ClipboardEvent('paste', { clipboardData: data, bubbles: true, cancelable: true }));
    JS

    assert_equal %w[AAAACCCC GGGGTTTT], all('#samples input[data-col="sequence"]').map(&:value)
  end

  test 'a saved oligo can be marked for removal, restored, and removed on save' do
    primer_set = edit_set_with(%w[keep_1_LEFT drop_1_LEFT])
    row = find('#samples tr.nested-fields') { |tr| tr.find('input[data-col="name"]').value == 'drop_1_LEFT' }
    row.find('[data-nested-remove]').click

    assert row.matches_css?('.is-removed')

    row.find('[data-nested-remove]').click

    assert_not row.matches_css?('.is-removed')

    row.find('[data-nested-remove]').click
    click_on 'Update Primer set'

    assert_text 'successfully updated'
    assert_equal %w[keep_1_LEFT], primer_set.oligos.reload.map(&:name)
  end

  test 'a FASTA upload on the edit page updates oligos by name and adds the rest' do
    fasta = Rails.root.join('tmp/system_test_merge.fasta')
    File.write(fasta, ">old_1_LEFT\nTTTTTTTTTT\n>new_1_RIGHT\nGGGGGGGGGG\n")
    edit_set_with(%w[old_1_LEFT])
    attach_file 'fasta_upload', fasta, make_visible: true

    assert_selector '#samples .nested-fields', count: 2
    assert_equal({ 'old_1_LEFT' => 'TTTTTTTTTT', 'new_1_RIGHT' => 'GGGGGGGGGG' },
                 all('#samples tr.nested-fields').to_h do |tr|
                   [tr.find('input[data-col="name"]').value, tr.find('input[data-col="sequence"]').value]
                 end)
  ensure
    FileUtils.rm_f(fasta)
  end
end
