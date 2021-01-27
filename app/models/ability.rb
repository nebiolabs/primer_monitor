# frozen_string_literal: true

class Ability
  include CanCan::Ability

  # See the cancancan wiki for syntax details:
  # https://github.com/CanCanCommunity/cancancan/blob/develop/docs/Defining-Abilities.md

  def initialize(user)
    clear_aliased_actions
    alias_action :show, :show, to: :read
    alias_action :edit, to: :update

    @user = user || User.new

    guest_ability
    return if @user.new_record?

    admin_ability if @user.has_role?(:administrator)
  end

  def guest_ability
    can :show, WelcomeController
    can :show, AboutController
    can :show, HelpController
    can :show, HistoryController
    can :show, UserEmailConfirmationsController
    can :create, UserSession
    can %i[new create], User
    can :read, PrimerSet
    can :read, Oligo
  end

  def admin_ability
    can :manage, :all
  end
end
