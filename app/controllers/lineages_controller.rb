# frozen_string_literal: true

class LineagesController < ApplicationController

  def index
    authorize! :index, LineagesController
    if params.key? :organism_name
      @organism = Organism.find_by(name: params[:organism_name])
      @lineages = Lineage.where(organism_id: @organism.id)
    else
      redirect_to lineage_variants_url, status: 301
    end
  end

  def show
    authorize! :show, LineagesController
    @lineage = Lineage.find_by(name: params[:name])
  end

  def to_param
    name
  end
end
