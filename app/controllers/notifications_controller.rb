# frozen_string_literal: true

# Admin dashboard of the notification emails sent to (or skipped for) users
class NotificationsController < ApplicationController
  before_action { authorize! :manage, VerifiedNotification }

  def index
    set_filters
    period = @from.beginning_of_day..@to.end_of_day
    notifications = VerifiedNotification.filtered(period, status: @status, user_query: @query)
    @status_counts = notifications.group(:status).count
    @recipient_count = notifications.where(status: 'Sent').distinct.count(:user_id)
    @notifications = notifications.with_alert_counts.preload(:user).order(updated_at: :desc)
  end

  def show
    @notification = VerifiedNotification.includes(:user).find(params.expect(:id))
    @alerts = @notification.proposed_notifications
                           .includes(:primer_set, :oligo, :detailed_geo_location_alias)
                           .joins(:primer_set, :oligo)
                           .order('primer_sets.name, proposed_notifications.coordinate, oligos.name')
  end

  private

  # the last 30 days of all notifications unless the dashboard's filter form says otherwise
  def set_filters
    filters = params.permit(:from, :to, :status, :q)
    @from = parse_date(filters[:from]) || 30.days.ago.to_date
    @to = parse_date(filters[:to]) || Date.current
    @status = filters[:status].presence_in(VerifiedNotification::STATUSES)
    @query = filters[:q].to_s.strip
  end

  def parse_date(value)
    Date.iso8601(value) if value.present?
  rescue Date::Error
    nil
  end
end
