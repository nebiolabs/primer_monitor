# frozen_string_literal: true

class LineagesController < ApplicationController
  def index
    authorize! :index, LineagesController
  end
end
