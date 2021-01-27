# frozen_string_literal: true

# stores information about users of this system (including submitters and administrators/operators)
class User < ApplicationRecord
  acts_as_authentic do |c|
    c.crypto_provider = ::Authlogic::CryptoProviders::SCrypt
  end

  has_many :primer_sets
  has_many :user_roles, dependent: :destroy
  has_many :roles, through: :user_roles
  accepts_nested_attributes_for :user_roles, reject_if: :all_blank, allow_destroy: true

  validates :email,
            format: {
              with: /@/,
              message: 'should look like an email address.'
            },
            length: { maximum: 100 },
            uniqueness: {
              case_sensitive: false,
              if: :will_save_change_to_email?
            },
            allow_blank: true

  validates :login,
            length: { within: 3..100 },
            uniqueness: {
              case_sensitive: false,
              if: :will_save_change_to_login?
            },
            allow_blank: true

  validates :password,
            confirmation: { if: :require_password? },
            length: {
              minimum: 8,
              if: :require_password?
            }
  validates :password_confirmation,
            length: {
              minimum: 8,
              if: :require_password?
            }

  before_validation :set_login_from_email

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

  def has_role?(role_to_test)
    role_to_test_ary = if role_to_test.is_a?(Array)
                         role_to_test
                       else
                         [role_to_test]
                       end
    !(role_to_test_ary & role_symbols).empty?
  end

  def deliver_password_reset_instructions!
    reset_perishable_token!

    Notifier.deliver_password_reset_instructions(self)
  end

  def confirm!
    self.confirmed = true
    save
  end
end
