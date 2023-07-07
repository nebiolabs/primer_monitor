# frozen_string_literal: true

class LineagesController < ApplicationController

  def index
    authorize! :index, LineagesController
    @organism = Organism.where(ncbi_taxon_id: 2697049)
  end

  def show
    authorize! :show, LineagesController
  end
end
