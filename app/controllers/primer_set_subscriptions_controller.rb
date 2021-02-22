# frozen_string_literal: true

class PrimerSetSubscriptionsController < ApplicationController
  load_and_authorize_resource

  def create
    # if the user is subscribing, they must want email...
    current_user.send_primer_updates || current_user.update_column(:send_primer_updates, true)

    @primer_set_subscription = PrimerSetSubscription.find_or_initialize_by(primer_set_subscription_params)
    respond_to do |format|
      format.js if @primer_set_subscription.save
    end
  end

  def destroy
    primer_set_subscription = PrimerSetSubscription.find(params[:id])
    @primer_set_id = primer_set_subscription.primer_set_id
    primer_set_subscription.destroy
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
