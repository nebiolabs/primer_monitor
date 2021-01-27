class AboutController < ApplicationController
  def show
    authorize! :show, AboutController
  end
end
