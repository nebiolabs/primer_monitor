# frozen_string_literal: true

# This class is used to set up functionality that is common to all controllers
class ApplicationController < ActionController::Base
  check_authorization :unless => :devise_controller?
  before_action :configure_permitted_parameters, if: :devise_controller?
  before_action :update_last_request_time


  rescue_from CanCan::AccessDenied do |exception|
    respond_to do |format|
      format.json { head :forbidden, content_type: 'text/html' }
      if current_user
        format.html { redirect_to main_app.root_url, notice: exception.message }
      else
        format.html { redirect_to new_user_session_path, notice: "Please log in to gain access" }
      end
      format.js { head :forbidden, content_type: 'text/html' }
    end
  end

  protected

  def handle_unverified_request
    raise ActionController::InvalidAuthenticityToken
  end

  def configure_permitted_parameters
    devise_parameter_sanitizer.permit(:sign_up, keys: [:name])
    devise_parameter_sanitizer.permit(:account_update, keys: [:name])
  end

  def update_last_request_time
    current_user && current_user.update_attribute(:last_request_at, DateTime.now())
  end


end
