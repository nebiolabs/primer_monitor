# frozen_string_literal: true

class VerifiedNotification < ApplicationRecord
  belongs_to :user
  has_many :proposed_notifications

  # collects unsent notifications by user_id
  def self.find_or_create_verified_notifications!
    new_notifications = ProposedNotification.where(verified_notification_id: nil)
                                            .select(:user_id)
                                            .group(:user_id)

    new_vns = new_notifications.map do |notification|
      vn = VerifiedNotification.find_or_create_by!(user_id: notification.user_id, status: 'Unsent')
      vn.proposed_notifications += ProposedNotification.where(user_id: notification.user_id)
                                                       .where(verified_notification_id: nil)
      vn
    end
    existing_vns = VerifiedNotification.where.not(user_id: new_notifications.map(&:user_id))
                                       .where(status: 'Unsent')
    existing_vns + new_vns
  end
end
