# frozen_string_literal: true

# this is the landing page for the application
class PrimerVariantSummariesController < ApplicationController
  def show
    authorize! :show, PrimerVariantSummariesController
    @organism = Organism.find_by(slug: params[:organism_name])
  end
end
