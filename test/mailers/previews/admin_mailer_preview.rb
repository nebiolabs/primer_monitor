# frozen_string_literal: true

# Preview all emails at http://localhost:3000/rails/mailers/admin_mailer
class AdminMailerPreview < ActionMailer::Preview
  def notification_summary_email
    notifications = VerifiedNotification.where.not(status: 'Unsent').order(updated_at: :desc).limit(10).to_a
    AdminMailer.notification_summary_email(User.first, notifications, Date.current)
  end
end
