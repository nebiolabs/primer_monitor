# frozen_string_literal: true

class ProposedNotification < ApplicationRecord
  belongs_to :user
  belongs_to :oligo
  belongs_to :primer_set
  belongs_to :primer_set_subscription
  belongs_to :subscribed_geo_location
  belongs_to :verified_notification
  belongs_to :detailed_geo_location_alias

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

    IdentifyPrimersForNotification.includes(:detailed_geo_location).all.each do |primer_record|
      pn = construct_notification_record(primer_record)

      potential_notifications << pn unless existing_notification_cache.key?(pn.cache_key)
    end

    potential_notifications
  end

  def self.construct_notification_record(primer_record)
    subscribed_alias = JoinSubscribedLocationToId.find_by(user_id: primer_record.user_id,
                                                          detailed_geo_location_id:
                                                            primer_record.detailed_geo_location_id)
    subscribed_geo_location_id = SubscribedGeoLocation.find_by(user_id: primer_record.user_id,
                                                               detailed_geo_location_alias_id:
                                                                 subscribed_alias.detailed_geo_location_alias_id).id

    primer_set_subscription_id = PrimerSetSubscription.find_by(user_id: primer_record.user_id,
                                                               primer_set_id: primer_record.primer_set_id).id

    ProposedNotification.new(primer_set_id: primer_record.primer_set_id, user_id: primer_record.user_id,
                             oligo_id: primer_record.oligo_id, coordinate: primer_record.coords,
                             subscribed_geo_location_id: subscribed_geo_location_id,
                             primer_set_subscription_id: primer_set_subscription_id,
                             detailed_geo_location_alias_id: primer_record.detailed_geo_location
                                                                          .detailed_geo_location_alias_id,
                             fraction_variant: primer_record.fraction_variant)
  end
end
