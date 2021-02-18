# frozen_string_literal: true

class PrimerSet < ApplicationRecord
  belongs_to :user
  belongs_to :organism
  has_many :oligos, dependent: :destroy

  enum status: { pending: 'pending', complete: 'complete', failed: 'failed' }

  accepts_nested_attributes_for :oligos, reject_if: :all_blank, allow_destroy: true

end
