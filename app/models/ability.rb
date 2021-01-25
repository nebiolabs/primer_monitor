# frozen_string_literal: true

class Ability
  include CanCan::Ability

  # See the cancancan wiki for syntax details:
  # https://github.com/CanCanCommunity/cancancan/blob/develop/docs/Defining-Abilities.md

  def initialize(user)
    clear_aliased_actions
    alias_action :index, :show, to: :read
    alias_action :edit, to: :update

    @user = user || User.new

    guest_ability
    return if @user.new_record?

    admin_ability if @user.has_role?(:administrator)
  end

  def guest_ability
    can :index, WelcomeController
    can :create, UserSession
  end

  def admin_ability
    can :manage, :all
  end
end
