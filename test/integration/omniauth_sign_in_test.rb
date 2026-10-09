# frozen_string_literal: true

require 'test_helper'

class OmniauthSignInTest < ActionDispatch::IntegrationTest
  setup do
    # Devise gives OmniAuth its /users/auth path prefix when routes are drawn, and Rails loads routes lazily in test,
    # so a POST that is a process's first request would otherwise miss the OmniAuth middleware.
    Rails.application.reload_routes_unless_loaded
    OmniAuth.config.test_mode = true
  end

  teardown do
    OmniAuth.config.test_mode = false
    OmniAuth.config.mock_auth.except!(:google_oauth2, :entra_id)
  end

  test 'the log in page offers each provider as a POST button that bypasses turbo' do
    get new_user_session_path

    %w[google_oauth2 entra_id].each do |provider|
      assert_select "form[method=post][action='/users/auth/#{provider}'][data-turbo=false] button"
    end
  end

  %w[google_oauth2 entra_id].each do |provider|
    test "#{provider} sign-in signs in the account with that email" do
      OmniAuth.config.mock_auth[provider.to_sym] =
        OmniAuth::AuthHash.new(provider:, uid: 'uid-1', info: { email: users(:one).email, name: 'Frank Gehry' })

      post "/users/auth/#{provider}"
      follow_redirect! # to the callback

      assert_redirected_to root_path
      assert_equal users(:one).id, session['warden.user.user.key']&.first&.first
    end
  end

  test 'an existing user can still log in with their password after signing in with Microsoft' do
    users(:one).update!(password: 'local-password-123')
    OmniAuth.config.mock_auth[:entra_id] =
      OmniAuth::AuthHash.new(provider: 'entra_id', uid: 'uid-1', info: { email: users(:one).email })
    post '/users/auth/entra_id'
    follow_redirect!
    delete destroy_user_session_path

    post user_session_path, params: { user: { email: users(:one).email, password: 'local-password-123' } }

    assert_equal users(:one).id, session['warden.user.user.key']&.first&.first
  end

  test 'a failed sign-in returns to the log in page' do
    OmniAuth.config.mock_auth[:entra_id] = :invalid_credentials

    post '/users/auth/entra_id'
    follow_redirect! # to the failure endpoint

    assert_redirected_to new_user_session_path
  end
end
