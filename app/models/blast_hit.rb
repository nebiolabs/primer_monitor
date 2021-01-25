# frozen_string_literal: true

class BlastHit < ApplicationRecord
  belongs_to :organism
  belongs_to :oligo
end
