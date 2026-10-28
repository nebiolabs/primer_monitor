# frozen_string_literal: true

require 'test_helper'

class LineagesControllerTest < ActionDispatch::IntegrationTest
  setup do
    @organism = organisms(:sars_cov2)
    taxon = organism_taxa(:sars_cov_2_taxon)
    primer_sets(:one).update!(status: 'complete')
    @lineage = Lineage.create!(name: 'XFG.CTRL', organism: @organism)
    call = LineageCall.create!(taxon: 'XFG_CTRL_ACC', lineage: @lineage)
    12.times do |i|
      record = FastaRecord.create!(strain: "ctrl/#{i}", organism_taxon_id: taxon.id, lineage_call: call,
                                   detailed_geo_location: detailed_geo_locations(:darwin), date_collected: Date.current)
      VariantSite.create!(fasta_record_id: record.id, organism_taxon_id: taxon.id, ref: nil,
                          ref_start: 110, ref_end: 111, variant_type: 'X', variant: 'T')
    end
    ActiveRecord::Base.connection.execute('REFRESH MATERIALIZED VIEW lineage_info')
  end

  test 'the lineage page shows where the lineage is seen and its primer sets' do
    get organism_lineage_url(organism_slug: @organism.slug, name: @lineage.name)

    assert_response :success
    assert_select 'div#lineage_map[data-points*="Australia"]'
    assert_select 'tr.is-clickable[data-map-point="Australia"]'
    assert_select 'p#lineage_config', text: /primer_set_tracks/
    assert_select 'td', text: primer_sets(:one).name
  end

  test 'serves a primer set as a BED track' do
    get organism_lineage_primer_track_url(organism_slug: @organism.slug, lineage: @lineage.name,
                                          set: primer_sets(:one).id, format: :text)

    assert_response :success
    assert_includes response.body, '95'
  end

  test 'returns 404 for an unknown primer set' do
    get organism_lineage_primer_track_url(organism_slug: @organism.slug, lineage: @lineage.name, set: 0,
                                          format: :text)

    assert_response :not_found
  end
end
