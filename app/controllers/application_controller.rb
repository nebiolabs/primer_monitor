# frozen_string_literal: true

# This class is used to set up functionality that is common to all controllers
class ApplicationController < ActionController::Base
  check_authorization

  helper_method :current_user_session, :current_user

  rescue_from CanCan::AccessDenied do |exception|
    respond_to do |format|
      format.json { head :forbidden, content_type: 'text/html' }
      if current_user
        format.html { redirect_to main_app.root_url, notice: exception.message }
      else
        format.html { redirect_to new_user_sessions_path, notice: "Please log in to gain access" }
      end
      format.js { head :forbidden, content_type: 'text/html' }
    end
  end

  private

  def current_user_session
    return @current_user_session if defined?(@current_user_session)

    @current_user_session = UserSession.find
  end

  def current_user
    return @current_user if defined?(@current_user)

    @current_user = current_user_session&.user
  end
end
