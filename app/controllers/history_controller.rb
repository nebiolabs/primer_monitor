# frozen_string_literal: true

class HistoryController < ApplicationController
  def show
    authorize! :show, HistoryController
  end
end
