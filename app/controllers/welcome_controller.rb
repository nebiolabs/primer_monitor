# frozen_string_literal: true

# this is the landing page for the application
class WelcomeController < ApplicationController
  def index
    authorize! :index, WelcomeController
  end
end
