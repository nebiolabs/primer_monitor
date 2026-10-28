# frozen_string_literal: true

require 'test_helper'

module Lineages
  class GeographyTest < ActiveSupport::TestCase
    setup do
      @organism = organisms(:sars_cov2)
      @taxon = organism_taxa(:sars_cov_2_taxon)
      @lineage = Lineage.create!(name: 'XFG.GEO', organism: @organism)
      @call = LineageCall.create!(taxon: 'XFG_GEO_ACC', lineage: @lineage)
    end

    test 'counts the lineage share of each map position in the window' do
      FastaRecord.create!(strain: 'geo/in', organism_taxon_id: @taxon.id, lineage_call: @call,
                          detailed_geo_location: detailed_geo_locations(:darwin), date_collected: Date.current)
      FastaRecord.create!(strain: 'geo/other', organism_taxon_id: @taxon.id,
                          detailed_geo_location: detailed_geo_locations(:darwin), date_collected: Date.current)
      FastaRecord.create!(strain: 'geo/old', organism_taxon_id: @taxon.id, lineage_call: @call,
                          detailed_geo_location: detailed_geo_locations(:darwin), date_collected: 200.days.ago.to_date)

      result = Geography.call(organism: @organism, lineage: @lineage, since: Date.current - 30)

      assert_equal 1, result[:lineage_sequences]
      assert_equal 2, result[:sequences]
      point = result[:points].sole

      assert_equal 'Australia', point[:name]
      assert_in_delta 0.5, point[:fraction]
    end
  end
end
