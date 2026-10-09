# frozen_string_literal: true

class PrimerSetSubscription < ApplicationRecord
  belongs_to :user
  belongs_to :primer_set
  has_many :proposed_notifications, dependent: :destroy

  validates :user_id, uniqueness: { scope: :primer_set_id }

  Summary = Struct.new(:subscription, :overlapping_variants, :current_alerts, :notifications_sent, :last_notified_at,
                       keyword_init: true) do
    delegate :primer_set, to: :subscription
  end

  # The user's active subscriptions, each with how many distinct variants overlapped the set's primers within the
  # user's lookback window, how many overlaps currently pass the user's alert settings, and what was emailed.
  def self.summaries_for(user)
    subscriptions = where(user:, active: true).joins(:primer_set).includes(primer_set: :organism)
                                              .order('primer_sets.name')
    ids = subscriptions.map(&:primer_set_id)
    overlaps = overlapping_variant_counts(ids, user.lookback_days)
    alerts = current_alert_counts(user, ids)
    sent = ProposedNotification.sent.where(user:, primer_set_id: ids).group(:primer_set_id)
    sent_counts = sent.count
    last_sent = sent.maximum('verified_notifications.updated_at')

    subscriptions.map do |subscription|
      id = subscription.primer_set_id
      Summary.new(subscription:, overlapping_variants: overlaps[id], current_alerts: alerts[id],
                  notifications_sent: sent_counts.fetch(id, 0), last_notified_at: last_sent[id])
    end
  end

  # primer_set_id => distinct variants overlapping its primers in the last lookback_days ({} until the
  # backend has populated the overlap views, e.g. in development)
  def self.overlapping_variant_counts(primer_set_ids, lookback_days)
    return {} if primer_set_ids.empty? || !materialized_view_populated?('oligo_variant_overlaps')

    connection.select_rows(sanitize_sql([<<~SQL.squish, primer_set_ids, lookback_days])).to_h
      SELECT primer_set_id, count(DISTINCT variant_id) FROM oligo_variant_overlaps
      WHERE primer_set_id IN (?) AND date_collected >= current_date - ?::integer
      GROUP BY primer_set_id
    SQL
  end

  def self.current_alert_counts(user, primer_set_ids)
    return {} if primer_set_ids.empty? || !materialized_view_populated?('identify_primers_for_notifications')

    IdentifyPrimersForNotification.where(user_id: user.id, primer_set_id: primer_set_ids).group(:primer_set_id).count
  end

  def self.materialized_view_populated?(name)
    connection.select_value(sanitize_sql(['SELECT ispopulated FROM pg_matviews WHERE matviewname = ?', name]))
  end

  # creates a hash of primer_set_ids -> PrimerSetSubscriptions for the specified user
  def self.subscriptions_for_user_by_primer_set(user)
    return {} unless user

    PrimerSetSubscription.where(user_id: user.id, active: true).index_by(&:primer_set_id)
  end
end
