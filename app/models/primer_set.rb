# frozen_string_literal: true

class PrimerSet < ApplicationRecord
  belongs_to :user
  has_many :oligos
end
