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
    # Either create a User record or update it based on the provider (Google) and the UID
    where(provider: auth.provider, uid: auth.uid).first_or_create do |user|
      user.token = auth.credentials.token
      user.expires = auth.credentials.expires
      user.expires_at = auth.credentials.expires_at
      user.refresh_token = auth.credentials.refresh_token
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
