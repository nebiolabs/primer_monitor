# frozen_string_literal: true

class LineageCaller < ApplicationRecord
  has_many :organism_taxa, dependent: :nullify
end
