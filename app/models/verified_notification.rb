# frozen_string_literal: true

class VerifiedNotification < ApplicationRecord
  belongs_to :user
  has_many :proposed_notifications, dependent: :destroy

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

  # Records every pending notification as handled without emailing it, except for the given addresses.
  # Run after changes that would otherwise flood subscribers with a backlog (see notifications:baseline).
  def self.skip_unsent!(except_emails: [])
    keep = except_emails.map { |email| email.strip.downcase }
    find_or_create_verified_notifications!.reject { |vn| keep.include?(vn.user.email.downcase) }
                                          .each { |vn| vn.update!(status: 'Skipped') }
  end
end
