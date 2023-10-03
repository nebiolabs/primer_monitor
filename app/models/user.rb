# frozen_string_literal: true

# stores information about users of this system (including submitters and administrators)
class User < ApplicationRecord
  devise :database_authenticatable, :registerable, :recoverable, :rememberable,
         :validatable, :confirmable, :omniauthable, omniauth_providers: [:google_oauth2]

  has_many :user_roles, dependent: :destroy
  has_many :roles, through: :user_roles
  has_many :subscribed_geo_locations, dependent: :destroy, inverse_of: :user
  has_many :detailed_geo_location_aliases, through: :subscribed_geo_locations
  has_many :primer_set_subscriptions, dependent: :destroy
  has_many :primer_sets, through: :primer_set_subscriptions
  has_many :verified_notifications, dependent: :destroy

  accepts_nested_attributes_for :user_roles, reject_if: :all_blank, allow_destroy: true

  before_validation :set_login_from_email

  def self.from_omniauth(auth)
    data = auth.info
    Rails.logger.debug("Attempting to log in via oauth with data: #{data}")
    user_attribs = {
      email: data['email'], first: data['first_name'],
      last: data['last_name'],
      password: Devise.friendly_token[0, 20]
    }
    User.create_with(user_attribs).find_or_create_by!(email: data['email']) do |user|
      Rails.logger.info("Creating new user using : #{user_attribs}")
      user.skip_confirmation!
    end
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
