# frozen_string_literal: true

class LineagesController < ApplicationController
  def index
    authorize! :index, LineagesController
    if params.key? :organism_name
      @organism = Organism.find_by(slug: params[:organism_name])

      @lineages = LineageInfo.where(organism_id: @organism.id).select(:name, :times_seen, :last_seen).to_a
    else
      # hardcoded legacy redirect
      redirect_to organism_lineage_variants_url(Organism.find_by(slug: 'sars-cov-2').slug), status: :moved_permanently
    end
  end

  def show
    authorize! :show, LineagesController
    @organism = Organism.find_by(name: params[:organism_name])
    @lineage = Lineage.find_by(name: params[:name])
    @lineage_info = LineageInfo.find_by(organism_id: @organism.id, name: @lineage.name)
  end

  def to_param
    name
  end
end
