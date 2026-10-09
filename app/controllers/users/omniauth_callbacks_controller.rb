# frozen_string_literal: true

module Users
  # Sign-in callbacks from the identity providers configured in config/initializers/devise.rb
  class OmniauthCallbacksController < Devise::OmniauthCallbacksController
    PROVIDER_NAMES = { google_oauth2: 'Google', entra_id: 'Microsoft' }.freeze

    def google_oauth2
      sign_in_from_omniauth
    end

    def entra_id
      sign_in_from_omniauth
    end

    protected

    def after_omniauth_failure_path_for(_scope)
      new_user_session_path
    end

    def after_sign_in_path_for(resource_or_scope)
      stored_location_for(resource_or_scope) || root_path
    end

    private

    def sign_in_from_omniauth
      auth = request.env['omniauth.auth']
      @user = User.from_omniauth(auth)
      sign_in_and_redirect @user, event: :authentication
      set_flash_message(:notice, :success, kind: PROVIDER_NAMES.fetch(auth.provider.to_sym)) if is_navigational_format?
    rescue ActiveRecord::RecordInvalid => e
      Rails.logger.warn("#{auth&.provider} sign-in failed: #{e.message}")
      redirect_to new_user_session_path, alert: "Could not sign in: #{e.record.errors.full_messages.to_sentence}"
    end
  end
end
