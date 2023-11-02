# frozen_string_literal: true

class OrganismTaxon < ApplicationRecord
  belongs_to :lineage_caller
end
