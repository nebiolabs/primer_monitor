# frozen_string_literal: true

class PrimerSet < ApplicationRecord
  belongs_to :user
  belongs_to :organism
  belongs_to :amplification_method
  has_many :oligos, -> { order('oligos.ref_start') }, dependent: :destroy, inverse_of: :primer_set
  has_many :subscriptions, dependent: :destroy, class_name: 'PrimerSetSubscription'
  has_many :subscribers, through: :subscriptions, source: :user

  enum status: { pending: 'pending', complete: 'complete', failed: 'failed' }

  accepts_nested_attributes_for :oligos, reject_if: :all_blank, allow_destroy: true

  def to_s
    name
  end

  def display_url
    citation_url || doi_url
  end

  def doi_url
    "https://doi.org/#{doi}" unless doi.blank?
  end

  def subscription_for_user(user)
    subscriptions.where(user_id: user.id).first
  end
end
