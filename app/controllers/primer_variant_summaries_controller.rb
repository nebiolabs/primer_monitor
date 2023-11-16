# frozen_string_literal: true

# this is the landing page for the application
class PrimerVariantSummariesController < ApplicationController
  def index
    @organism = Organism.find_by(slug: params[:organism_name])
    authorize! :show, PrimerVariantSummariesController
  end
end
