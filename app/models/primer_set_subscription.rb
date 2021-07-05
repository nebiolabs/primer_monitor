# frozen_string_literal: true

class PrimerSetSubscription < ApplicationRecord
  belongs_to :user
  belongs_to :primer_set
  has_many :proposed_notifications, dependent: :destroy

  validates_uniqueness_of :user_id, scope: :primer_set_id

  # creates a hash of primer_set_ids -> PrimerSetSubscriptions for the specified user
  def self.subscriptions_for_user_by_primer_set(user)
    return {} unless user

    PrimerSetSubscription.where(user_id: user.id, active: true).each_with_object({}) do |pss, h|
      h[pss.primer_set_id] = pss
    end
  end
end
