# frozen_string_literal: true

class PrimerSet < ApplicationRecord
  belongs_to :user
  belongs_to :organism
  has_many :oligos, dependent: :destroy
  has_many :subscriptions, dependent: :destroy, class_name: 'PrimerSetSubscription'
  has_many :subscribers, through: :subscriptions, source: :user

  enum status: { pending: 'pending', complete: 'complete', failed: 'failed' }

  accepts_nested_attributes_for :oligos, reject_if: :all_blank, allow_destroy: true

  def subscribed?(user)
    return {} unless user

    { id => subscriptions.where(user_id: user.id) }
  end
end
