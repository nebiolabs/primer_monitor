# frozen_string_literal: true

# stores information about users of this system (including submitters and administrators)
class User < ApplicationRecord
  # omniauth providers are those configured in config/initializers/devise.rb
  devise :database_authenticatable, :registerable, :recoverable, :rememberable,
         :validatable, :confirmable, :omniauthable

  has_many :user_roles, dependent: :destroy
  has_many :roles, through: :user_roles
  has_many :subscribed_geo_locations, dependent: :destroy, inverse_of: :user
  has_many :detailed_geo_location_aliases, through: :subscribed_geo_locations
  has_many :primer_set_subscriptions, dependent: :destroy
  has_many :primer_sets, through: :primer_set_subscriptions
  has_many :verified_notifications, dependent: :destroy

  accepts_nested_attributes_for :user_roles, reject_if: :all_blank, allow_destroy: true

  before_validation :set_login_from_email

  # The user signing in through Google or Entra ID: the account already linked to that identity, else the
  # existing account with their email (linked from now on), else a new account.
  # The identity provider has verified the email, so the account counts as confirmed.
  def self.from_omniauth(auth)
    # Entra ID sends no email claim for accounts without a mail attribute; the sign-in name (UPN) is then the email
    auth.info.email = auth.info.nickname if auth.info.email.blank? && auth.info.nickname.to_s.include?('@')
    user = find_by(provider: auth.provider, uid: auth.uid) || link_by_email(auth) || build_from_omniauth(auth)
    user.skip_confirmation! unless user.confirmed?
    user.save!
    user
  end

  def self.link_by_email(auth)
    email = auth.info.email.to_s.strip.downcase
    return if email.empty?

    find_by('lower(email) = ?', email)&.tap { |user| user.assign_attributes(provider: auth.provider, uid: auth.uid) }
  end

  def self.build_from_omniauth(auth)
    info = auth.info
    Rails.logger.info("Creating new #{auth.provider} user for #{info.email}")
    first, last = info.name.to_s.split(' ', 2)
    new(provider: auth.provider, uid: auth.uid, email: info.email,
        first: info.first_name.presence || first.presence || info.email,
        last: info.last_name.presence || last.presence || '-',
        password: Devise.friendly_token[0, 20])
  end

  def subscribed_detailed_geo_location_alias_ids
    subscribed_geo_locations.map(&:detailed_geo_location_alias_id)
  end

  def subscribed_detailed_geo_location_alias_ids=(dga_ids)
    recs = []
    dga_ids.each do |dga_id|
      next if dga_id.blank?

      recs << SubscribedGeoLocation.new(user_id: id, detailed_geo_location_alias_id: dga_id)
    end
    self.subscribed_geo_locations = recs
  end

  def to_s
    if first && last
      "#{first} #{last}"
    else
      login
    end
  end

  def set_login_from_email
    self.login ||= email
  end

  def role_symbols
    roles.map { |r| r.name.gsub(/\s+/, '_').downcase.to_sym }
  end

  def role?(role_to_test)
    role_to_test_ary = if role_to_test.is_a?(Array)
                         role_to_test
                       else
                         [role_to_test]
                       end
    !(role_to_test_ary & role_symbols).empty?
  end

  def formatted_email
    m = Mail::Address.new email
    m.display_name = "#{first.capitalize} #{last.capitalize}"
    m.format
  end
end
