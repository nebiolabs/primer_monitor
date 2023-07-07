# frozen_string_literal: true

class LineagesController < ApplicationController

  def index
    authorize! :index, LineagesController
    @organism = Organism.find_by(ncbi_taxon_id: 2697049)
    @lineages = Lineage.where(organism_id: @organism.id)
  end

  def show
    authorize! :show, LineagesController
  end
end
