# frozen_string_literal: true

class PrimerSetSubscriptionsController < ApplicationController
  load_and_authorize_resource

  def create
    primer_set = PrimerSet.find(params.expect(:primer_set_id)) # a stale or tampered id is a 404
    PrimerSetSubscription.find_or_initialize_by(user: current_user, primer_set:).update!(active: true)
    current_user.subscribe_to_primer_updates!
    redirect_back_or_to primer_set, notice: 'Subscribed.'
  end

  def destroy
    primer_set_subscription = PrimerSetSubscription.find(params.expect(:id))
    # setting this boolean to false is always going to be fine
    # rubocop:disable-next Rails/SkipsModelValidations
    primer_set_subscription.update_column(:active, false)
    redirect_back_or_to edit_user_registration_path, notice: 'Unsubscribed.'
  end

  # Only allow a list of trusted parameters through.
  def primer_set_subscription_params
    # to avoid non-user generated subscriptions, this always adds current_user
    params.permit(:id, :primer_set_id).merge(user_id: current_user.id)
  end
end
