# frozen_string_literal: true

require 'application_system_test_case'

# The new primer set form: adding and removing oligo rows (nested_fields.js) and filling them from a FASTA file.
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
end
