# frozen_string_literal: true

class AdminMailer < ApplicationMailer
  default from: "NEB Primer Monitor <#{ENV.fetch('ADMIN_EMAIL', nil)}>"

  # what the daily notifications run emailed (or skipped), linking to the dashboard for that day
  def notification_summary_email(admin, notifications, date)
    @date = date
    @sent = notifications.select { |vn| vn.status == 'Sent' }
    @skipped = notifications.select { |vn| vn.status == 'Skipped' }
    @alert_counts = ProposedNotification.where(verified_notification_id: notifications.map(&:id))
                                        .group(:verified_notification_id).count
    @url = notifications_url(from: date, to: date)
    mail(to: admin.formatted_email,
         subject: "Primer Monitor: #{@sent.size} notification #{'email'.pluralize(@sent.size)} sent #{date}")
  end
end
