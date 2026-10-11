# frozen_string_literal: true

class VerifiedNotification < ApplicationRecord
  belongs_to :user
  has_many :proposed_notifications, dependent: :destroy

  STATUSES = %w[Sent Skipped Unsent].freeze

  # each notification email with the number of alerts (proposed notifications) it carries
  scope :with_alert_counts, lambda {
    left_joins(:proposed_notifications).group(:id).select('verified_notifications.*',
                                                          'count(proposed_notifications.id) AS alert_count')
  }

  # last changed (sent, skipped or queued) within the given time range, optionally with the given status and for
  # users whose name or email contains user_query
  def self.filtered(range, status: nil, user_query: nil)
    notifications = where(updated_at: range)
    notifications = notifications.where(status:) if status
    return notifications if user_query.blank?

    notifications.joins(:user).where('users.email ILIKE :q OR users.first ILIKE :q OR users.last ILIKE :q',
                                     q: "%#{sanitize_sql_like(user_query)}%")
  end

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

  # emails the user their pending notifications and records the delivery; inactive users get no email, and their
  # notifications are recorded as skipped so reactivating them doesn't send the backlog
  def deliver!
    return update!(status: 'Skipped') unless user.active?

    Rails.logger.info("Sending notification to #{user.formatted_email}")
    notifications = ProposedNotification.where(verified_notification_id: id)
                                        .includes(:oligo, :primer_set,
                                                  :subscribed_geo_location, :detailed_geo_location_alias)
                                        .reorder('primer_sets.name, oligos.locus, oligos.category,' \
                                                 'region, subregion, division, subdivision')
    PrimerSetMailer.primer_overlap_notification_email(user_id, notifications).deliver_now
    update!(status: 'Sent')
  end

  # Records every pending notification as handled without emailing it, except for the given addresses.
  # Run after changes that would otherwise flood subscribers with a backlog (see notifications:baseline).
  def self.skip_unsent!(except_emails: [])
    keep = except_emails.map { |email| email.strip.downcase }
    find_or_create_verified_notifications!.reject { |vn| keep.include?(vn.user.email.downcase) }
                                          .each { |vn| vn.update!(status: 'Skipped') }
  end

  # emails each active administrator a summary of the given notifications, if there are any
  def self.send_summary!(notifications, date: Date.current)
    handled = notifications.reject { |vn| vn.status == 'Unsent' }
    return if handled.empty?

    User.joins(:roles).where(roles: { name: 'administrator' }, active: true).find_each do |admin|
      AdminMailer.notification_summary_email(admin, handled, date).deliver_now
    end
  end
end
