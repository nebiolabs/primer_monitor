# frozen_string_literal: true

class Organism < ApplicationRecord
  has_many :blast_hits
  has_many :amplicons
end
