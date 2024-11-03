# frozen_string_literal: true

class OrganismTaxon < ApplicationRecord
  belongs_to :lineage_caller
  belongs_to :organism

  def url
    "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=#{ncbi_taxon_id}"
  end
end
