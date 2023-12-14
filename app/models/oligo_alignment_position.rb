# frozen_string_literal: true

class OligoAlignmentPosition < ApplicationRecord
  belongs_to :oligo
  belongs_to :organism_taxon
end
