# frozen_string_literal: true

class PrimerSetMailer < ApplicationMailer
  default from: "NEB Primer Monitor <#{ENV['ADMIN_EMAIL']}>"

  def updated_primer_set_email(email, primer_set)
    @url = primer_set_url(primer_set)
    mail(to: email, subject: 'Primer Set Created or Updated')
  end

  def primer_overlap_notification_email(user_id, proposed_notifications)
    @user = User.find(user_id)
    @proposed_notifications = proposed_notifications
    mail(to: @user.formatted_email,
         subject: "#{@proposed_notifications.size} variants observed in your subscribed primer sets.")
  end
end
