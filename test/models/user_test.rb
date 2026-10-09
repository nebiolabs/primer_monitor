# frozen_string_literal: true

require 'test_helper'

class UserTest < ActiveSupport::TestCase
  test 'instantiation' do
    assert_not_nil User.new
  end

  test 'a new account is subscribed to World unless it chose locations' do
    defaulted = User.create!(first: 'A', last: 'B', email: 'a@example.org', password: 'password-123')
    chosen = User.create!(first: 'C', last: 'D', email: 'c@example.org', password: 'password-123',
                          subscribed_detailed_geo_location_alias_ids: [detailed_geo_location_aliases(:darwin).id])

    assert_equal [detailed_geo_location_aliases(:world)], defaulted.reload.detailed_geo_location_aliases
    assert_equal [detailed_geo_location_aliases(:darwin)], chosen.reload.detailed_geo_location_aliases
  end

  # --- User.from_omniauth (Google and Entra ID sign-in) ---

  def auth(email, provider: 'entra_id', uid: 'uid-1', **info)
    OmniAuth::AuthHash.new(provider:, uid:, info: { email:, name: 'Someone Else', **info })
  end

  test 'sign-in links an existing account with the same email, ignoring case' do
    user = User.from_omniauth(auth('Frank@Contemporary.example.org'))

    assert_equal users(:one), user
    assert_equal %w[entra_id uid-1], [user.reload.provider, user.uid]
    assert_equal user, User.from_omniauth(auth('someone-else@example.org')), 'later sign-ins find it by identity'
  end

  test 'sign-in uses the Entra sign-in name when there is no email claim' do
    user = User.from_omniauth(auth(nil, nickname: 'frank@contemporary.example.org'))

    assert_equal users(:one), user
  end

  test 'sign-in confirms an existing unconfirmed account' do
    users(:one).update_columns(confirmed_at: nil) # rubocop:disable Rails/SkipsModelValidations

    assert_predicate User.from_omniauth(auth('frank@contemporary.example.org')), :confirmed?
  end

  test 'sign-in creates a confirmed account for someone new' do
    user = assert_difference('User.count') do
      User.from_omniauth(auth('new@example.org', provider: 'google_oauth2', first_name: 'New', last_name: 'Person'))
    end

    assert_equal %w[New Person google_oauth2 uid-1], [user.first, user.last, user.provider, user.uid]
    assert_predicate user, :confirmed?
  end

  test 'sign-in falls back to the display name when first and last names are missing' do
    user = User.from_omniauth(auth('mononym@example.org', name: 'Cher'))

    assert_equal %w[Cher -], [user.first, user.last]
  end
end
