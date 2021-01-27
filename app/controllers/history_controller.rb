class HistoryController < ApplicationController
  def show
    authorize! :show, HistoryController
  end
end
