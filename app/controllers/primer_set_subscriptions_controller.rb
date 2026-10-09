# frozen_string_literal: true

class PrimerSetSubscriptionsController < ApplicationController
  load_and_authorize_resource

  def create
    current_user.subscribe_to_primer_updates!
    @primer_set_subscription = PrimerSetSubscription.find_or_initialize_by(primer_set_subscription_params)
    @primer_set_subscription.active = true
    respond_to do |format|
      format.js if @primer_set_subscription.save
      format.html do
        @primer_set_subscription.save!
        redirect_back_or_to @primer_set_subscription.primer_set, notice: 'Subscribed.'
      end
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
      format.html { redirect_back_or_to edit_user_registration_path, notice: 'Unsubscribed.' }
    end
  end

  # Only allow a list of trusted parameters through.
  def primer_set_subscription_params
    # to avoid non-user generated subscriptions, this always adds current_user
    params.permit(:id, :primer_set_id).merge(user_id: current_user.id)
  end
end
