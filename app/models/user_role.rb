# frozen_string_literal: true

# stores which users have which roles
class UserRole < ApplicationRecord
  belongs_to :user
  belongs_to :role

  validates :user, presence: true
  validates :role, presence: true
  validates :user_id, uniqueness: { scope: :role_id }
end
