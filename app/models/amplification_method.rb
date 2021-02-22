# frozen_string_literal: true

class AmplificationMethod < ApplicationRecord
  has_many :primer_sets

  def to_s
    name
  end
end
