# frozen_string_literal: true

class ProposedNotification < ApplicationRecord
  belongs_to :user
  belongs_to :oligo
  belongs_to :primer_set
  belongs_to :primer_set_subscription
  belongs_to :subscribed_geo_location

  UNIQUE_FIELDS = %i[primer_set_id user_id oligo_id coordinate
                     subscribed_geo_location_id primer_set_subscription_id].freeze

  def cache_key
    UNIQUE_FIELDS.map { |f| send f }.join
  end

  def self.existing_notification_cache
    @existing_notification_cache ||= ProposedNotification.pluck(:id, UNIQUE_FIELDS.join(','))
                                                         .each_with_object({}) do |pn_fields, h|
      h[pn_fields[1..].join] = pn_fields[0]
    end
  end

  def self.new_proposed_notifications
    potential_notifications = []

    IdentifyPrimersForNotification.all.each do |record|
      subscribed_alias = JoinSubscribedLocationToId.find_by(user_id: record.user_id,
                                                            detailed_geo_location_id:
                                                              record.detailed_geo_location_id)
      subscribed_geo_location_id = SubscribedGeoLocation.find_by(user_id: record.user_id,
                                                                 detailed_geo_location_alias_id:
                                                                   subscribed_alias.detailed_geo_location_alias_id).id

      detailed_geo_location_alias_id = DetailedGeoLocation.find_by(id: record.detailed_geo_location_id).detailed_geo_location_alias_id

      primer_set_subscription_id = PrimerSetSubscription.find_by(user_id: record.user_id,
                                                                 primer_set_id: record.primer_set_id).id

      pn = ProposedNotification.new(primer_set_id: record.primer_set_id, user_id: record.user_id,
                                    oligo_id: record.oligo_id, coordinate: record.coords,
                                    subscribed_geo_location_id: subscribed_geo_location_id,
                                    primer_set_subscription_id: primer_set_subscription_id,
                                    detailed_geo_location_alias_id: detailed_geo_location_alias_id,
                                    fraction_variant: record.fraction_variant)

      potential_notifications << pn unless existing_notification_cache.key?(pn.cache_key)
    end

    potential_notifications
  end

  def update_status!

  end
end
