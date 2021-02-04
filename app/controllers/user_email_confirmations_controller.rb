# frozen_string_literal: true

class UserEmailConfirmationsController < ApplicationController
  skip_authorization_check

  before_action :load_user_using_perishable_token

  def show
    if @user
      @user.confirm!
      flash[:notice] = 'Thank you for confirming your email address. Your account is now active.'
    end

    redirect_to root_url
  end

  private

  def load_user_using_perishable_token
    @user = User.find_using_perishable_token(params[:id], nil)
    flash[:notice] = 'Unable to confirm your account using the specified token.' unless @user
  end
end
