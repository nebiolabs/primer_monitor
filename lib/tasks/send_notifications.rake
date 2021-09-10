# frozen_string_literal: true

namespace :notifications do
  # sends messages and records successful delivery
  def send_notifications(verified_notifications)
    verified_notifications.each do |vn|
      Rails.logger.info("Sending notification to #{vn.user.formatted_email}")
      pns = ProposedNotification.where(verified_notification_id: vn.id)
                                .includes(:oligo, :primer_set,
                                          :subscribed_geo_location, :detailed_geo_location_alias)
                                .reorder('primer_sets.name, oligos.locus, oligos.category,' \
                                         'region, subregion, division, subdivision')
      PrimerSetMailer.primer_overlap_notification_email(vn.user_id, pns).deliver_now
      vn.status = 'Sent'
      vn.save!
    end
  end

  desc 'sends notifications about primer overlaps'
  task send: :environment do
    send_notifications(VerifiedNotification.find_or_create_verified_notifications!)
  end
end
