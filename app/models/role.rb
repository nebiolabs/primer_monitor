# frozen_string_literal: true

# represents the roles that users may take on
class Role < ApplicationRecord
  has_many :user_roles, dependent: :destroy
  has_many :users, through: :user_roles

  validates :name, uniqueness: true

  def to_s
    name
  end
end
