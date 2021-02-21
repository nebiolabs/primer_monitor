# frozen_string_literal: true

# stores information about users of this system (including submitters and administrators)
class User < ApplicationRecord
  devise :database_authenticatable, :registerable, :recoverable, :rememberable,
         :validatable, :confirmable, :omniauthable, omniauth_providers: [:google_oauth2]

  has_many :primer_sets
  has_many :user_roles, dependent: :destroy
  has_many :roles, through: :user_roles
  accepts_nested_attributes_for :user_roles, reject_if: :all_blank, allow_destroy: true

  before_validation :set_login_from_email

  def self.from_omniauth(auth)
    data = auth.info
    Rails.logger.debug("Attempting to log in via oauth with data: #{data}")
    User.find_or_create_by(email: data['email']) do |user|
      Rails.logger.info("Creating new user via oauth: #{data}")
      user.password = Devise.friendly_token[0, 20]
      user.first = user['first_name']
      user.last = user['last_name']
      user.confirmed = true
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
