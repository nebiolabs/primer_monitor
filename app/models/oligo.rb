# frozen_string_literal: true

class Oligo < ApplicationRecord
  belongs_to :primer_set
  has_many :blast_hits

  enum category: { F3: 'F3', FIP: 'FIP', LF: 'LF', LB: 'LB', BIP: 'BIP', B3: 'B3',
                   Forward: 'Forward', Probe: 'Probe', Reverse: 'Reverse' }

  validates :short_name, length: 1..5, allow_blank: true

  def to_s
    "#{name}: #{sequence}"
  end
end
