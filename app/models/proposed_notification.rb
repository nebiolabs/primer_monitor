# frozen_string_literal: true

class ProposedNotification < ApplicationRecord
  belongs_to :user
  belongs_to :oligo
  belongs_to :primer_set
  belongs_to :primer_set_subscription
  belongs_to :subscribed_geo_location
  
  def self.new_proposed_notifications()

    potential_notifications = []

    IdentifyPrimersForNotification.all.each do |record|

      oligo_id = Oligo.find_by(name: record.primer_name, primer_set_id: record.primer_set_id).id

      subscribed_alias = JoinSubscribedLocationToID.find_by(user_id: record.user_id, detailed_geo_location_id: record.detailed_geo_location_id)
      subscribed_geo_locations_id = SubscribedGeoLocation.find_by(user_id: record.user_id,
                                                                  detailed_geo_location_alias_id: subscribed_alias.detailed_geo_location_alias_id).id

      primer_set_subscriptions_id = PrimerSetSubscription.find_by(user_id: record.user_id, primer_set_id: record.primer_set_id).id

      return if ProposedNotification.exists?(primer_set_id: record.primer_set_id, user_id: record.user_id, oligo_id: oligo_id,
                                             coordinate: record.coords, subscribed_geo_locations_id: subscribed_geo_locations_id, 
                                             primer_set_subscriptions_id: primer_set_subscriptions_id)

      potential_notifications << ProposedNotification.new(primer_sets_id: record.primer_set_id, users_id: record.user_id, oligos_id: oligos_id,
                                                          coordinate: record.coords, subscribed_geo_locations_id: subscribed_geo_locations_id, 
                                                          primer_set_subscriptions_id: primer_set_subscriptions_id,
                                                          fraction_variant: record.fraction_variant)
                      
    end
    
    potential_notifications  
  end
end
  