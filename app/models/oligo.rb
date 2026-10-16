# frozen_string_literal: true

class Oligo < ApplicationRecord
  belongs_to :primer_set
  has_many :blast_hits, dependent: :destroy
  has_many :oligo_alignment_positions, dependent: :destroy
  has_many :lineage_variant_primer_overlaps, dependent: :destroy

  enum :category, { F3: 'F3', FIP: 'FIP', LF: 'LF', LB: 'LB', BIP: 'BIP', B3: 'B3',
                    Forward: 'Forward', Probe: 'Probe', Reverse: 'Reverse' }

  before_validation :infer_locus_and_category

  validates :short_name, length: 1..5, allow_blank: true

  validates :sequence, presence: true,
                       format: { with: /\A([ATUCGWSMKRYBDHVN]|\s)*\z/i, message: 'contains invalid bases' }

  validates :name, presence: true

  def to_s
    "#{name}: #{sequence}"
  end

  # name plus locus/category when set
  def label
    detail = [locus, category].compact_blank.join(' ')
    detail.present? ? "#{name} (#{detail})" : name
  end

  private

  # FASTA uploads carry only names; fill blanks from the name and leave anything entered by hand alone
  def infer_locus_and_category
    return if locus.present? && category.present?

    parsed = OligoNameParser.parse(name)
    self.locus = parsed[:locus] if locus.blank?
    self.category = parsed[:category] if category.blank?
  end
end
