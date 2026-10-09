# frozen_string_literal: true

require 'application_system_test_case'

# The igv.js genome viewer on the organism, primer set and lineage pages, against the fixture data server.
class GenomeViewerTest < ApplicationSystemTestCase
  # the lineage page reads lineage frequencies from this backend-refreshed view
  setup { ActiveRecord::Base.connection.execute('REFRESH MATERIALIZED VIEW lineage_info') }

  test 'the organism page shows the genes and variants tracks' do
    visit organism_path(organisms(:sars_cov2))

    assert_igv_tracks %w[Genes Variants]
  end

  test "a processed primer set's page shows its primers" do
    primer_sets(:two).update!(status: :complete)

    visit primer_set_path(primer_sets(:two))

    assert_igv_tracks %w[Genes CDC]
  end

  test 'changing the primer set selection on the lineage page updates the tracks' do
    visit organism_lineage_variants_path(organisms(:sars_cov2))

    assert_igv_tracks ['Genes', 'All Variants', 'CDC', 'Charité']

    page.execute_script(<<~JS)
      const picker = document.getElementById('primer_set_select').tomselect;
      picker.removeItem('Charite');
    JS

    assert_igv_tracks ['Genes', 'All Variants', 'CDC'], 'removing a primer set in the picker removes its track'
  end

  test 'moving between viewer pages with Turbo leaves one browser with that page’s tracks' do
    visit organism_path(organisms(:sars_cov2))

    assert_igv_tracks %w[Genes Variants]

    click_on 'Primer Status'

    assert_igv_tracks ['Genes', 'All Variants', 'CDC', 'Charité']

    go_back

    assert_igv_tracks %w[Genes Variants]
  end
end
