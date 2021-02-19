# frozen_string_literal: true

class AboutController < ApplicationController
  def show
    authorize! :show, AboutController
  end
end
