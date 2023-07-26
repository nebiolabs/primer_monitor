# frozen_string_literal: true

class PrimerSetSubscription < ApplicationRecord
  belongs_to :user
  belongs_to :primer_set
  has_many :proposed_notifications, dependent: :destroy

  validates :user_id, uniqueness: { scope: :primer_set_id }

  # creates a hash of primer_set_ids -> PrimerSetSubscriptions for the specified user
  def self.subscriptions_for_user_by_primer_set(user)
    return {} unless user

    PrimerSetSubscription.where(user_id: user.id, active: true).index_by(&:primer_set_id)
  end
end
