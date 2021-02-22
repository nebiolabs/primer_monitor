# frozen_string_literal: true

# stores information about users of this system (including submitters and administrators)
class User < ApplicationRecord
  devise :database_authenticatable, :registerable, :recoverable, :rememberable,
         :validatable, :confirmable, :omniauthable, omniauth_providers: [:google_oauth2]

  has_many :primer_sets
  has_many :user_roles, dependent: :destroy
  has_many :roles, through: :user_roles
  has_many :subscribed_geo_locations, dependent: :destroy
  has_many :geo_locations, through: :subscribed_geo_locations
  has_many :primer_set_subscriptions, dependent: :destroy
  has_many :primer_sets, through: :primer_set_subscriptions

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
    User.create_with(user_attribs).create_or_find_by!(email: data['email']) do |user|
      Rails.logger.info("Creating new user using : #{user_attribs}")
      user.skip_confirmation!
    end
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
end
