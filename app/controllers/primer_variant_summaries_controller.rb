# frozen_string_literal: true

# this is the landing page for the application
class PrimerVariantSummariesController < ApplicationController
  def index
    authorize! :index, PrimerVariantSummariesController
    @organism = Organism.find_by(slug: params[:organism_name])
  end
end
