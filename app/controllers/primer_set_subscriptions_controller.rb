# frozen_string_literal: true

class PrimerSetSubscriptionsController < ApplicationController
  load_and_authorize_resource

  def create
    # if the user is subscribing, they must want email...
    # setting this boolean to true is always going to be fine
    # rubocop:disable Rails/SkipsModelValidations
    current_user.send_primer_updates || current_user.update_column(:send_primer_updates, true)
    # rubocop:enable Rails/SkipsModelValidations
    @primer_set_subscription = PrimerSetSubscription.find_or_initialize_by(primer_set_subscription_params)
    @primer_set_subscription.active = true
    respond_to do |format|
      format.js if @primer_set_subscription.save
    end
  end

  def destroy
    primer_set_subscription = PrimerSetSubscription.find(params[:id])
    @primer_set_id = primer_set_subscription.primer_set_id
    # setting this boolean to false is always going to be fine
    # rubocop:disable Rails/SkipsModelValidations
    primer_set_subscription.update_column(:active, false)
    # rubocop:enable Rails/SkipsModelValidations
    respond_to do |format|
      format.js
    end
  end

  # Only allow a list of trusted parameters through.
  def primer_set_subscription_params
    # to avoid non-user generated subscriptions, this always adds current_user
    params.permit(:id, :primer_set_id).merge(user_id: current_user.id)
  end
end
