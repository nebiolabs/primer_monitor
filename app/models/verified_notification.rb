# frozen_string_literal: true

class VerifiedNotification < ApplicationRecord
    belongs_to :users
    has_many :proposed_notifications

    def self.group_notifications()

        # I want to find proposed notifications without a verified notification id, group by user, then
        # add a new verified record and put the new verified_record_id in the proposed notifications table for those
        # records
        # It seems to work, but there must be a better way of doing this
        new_notifications = ProposedNotification.where(:verified_notification_id => nil)
                                                .select(:user_id)
                                                .group(:user_id)
                                                
        new_notifications.each do |notification|

            unless VerifiedNotification.exists?(user_id: notification.user_id, status: 'Unsent')
                VerifiedNotification.new(user_id: notification.user_id, status: 'Unsent').save!
            end

            verified_id = VerifiedNotification.find_by(user_id: notification.user_id, status: 'Unsent').id

            ProposedNotification.where(:user_id => notification.user_id)
                                .where(:verified_notification_id => nil)
                                .update_all(:verified_notification_id => verified_id)
        end
    end
end