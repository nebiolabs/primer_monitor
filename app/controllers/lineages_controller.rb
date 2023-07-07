# frozen_string_literal: true

class LineagesController < ApplicationController

  def index
    authorize! :index, LineagesController
    @organism = Organism.find_by(ncbi_taxon_id: 2_697_049)
    @lineages = Lineage.where(organism_id: @organism.id)
  end

  def show
    authorize! :show, LineagesController
    @lineage = Lineage.find_by(id: params[:id])
  end
end
