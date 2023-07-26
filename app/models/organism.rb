# frozen_string_literal: true

class Organism < ApplicationRecord
  has_many :blast_hits, dependent: :destroy
  has_many :primer_sets, dependent: :destroy

  def to_s
    name
  end

  def full_name
    "#{name} (#{self.alias})"
  end
end
