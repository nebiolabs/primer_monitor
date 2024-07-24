# frozen_string_literal: true

class LineageVariantsController < ApplicationController
  def index
    authorize! :index, LineageVariantsController

    @organism = Organism.find_by(slug: params[:organism_name])

    @config = {
      "data_server": ENV['IGV_DATA_SERVER'],
      "organism_slug": @organism.slug,
      "organism_name": @organism.name
    }

    @track_data = @organism.lineage_variants_data(@config[:data_server], @config[:organism_slug])

    @default_lineage = params[:lineage] if params.key? 'lineage'

    @default_tracks = params[:primer_sets].split(',') if params.key? 'primer_sets'
  end
end
