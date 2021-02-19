# frozen_string_literal: true

class Organism < ApplicationRecord
  has_many :blast_hits
  has_many :primer_sets

  def to_s
    name
  end
end
