# frozen_string_literal: true

# provides login
class UserSessionsController < ApplicationController
  skip_authorization_check

  def new
    @user_session = UserSession.new
  end

  def create
    @user_session = UserSession.new(user_session_params.to_h)
    if @user_session.save
      redirect_to welcome_index_url
    else
      render action: :new
    end
  end

  def destroy
    current_user_session.destroy
    redirect_to new_user_sessions_url
  end

  private

  def user_session_params
    params.require(:user_session).permit(:login, :password, :remember_me)
  end
end
