# frozen_string_literal: true

class LineageVariantPrimerOverlap < ApplicationRecord
  belongs_to :organism
  belongs_to :oligo
end
