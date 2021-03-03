# frozen_string_literal: true

class VerifiedNotification < ApplicationRecord
    belongs_to :users
    has_many :proposed_notifications

    def group_notifications()

        # I want to find proposed notifications without a verified notification id, group by user, then
        # add a new verified record and put the new verified_record_id in the proposed notifications table for those
        # records
        new_notifications = ProposedNotification.find_by(verified_notifications_id: nil).group(:users_id)
                           
        # I think new_notifications is a hash? So this probably isn't the right syntax
        new_notifications.each do |notification|
            user_id = notifications.users_id
            VerifiedNotification.import(users_id: user_id, status: 'Unsent')
            notification.update(<new_verified_notification_id>)
        end

    end
end