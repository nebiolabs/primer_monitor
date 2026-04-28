# frozen_string_literal: true

class LineageVariantsController < ApplicationController
  def index
    authorize! :index, LineageVariantsController

    @organism = Organism.find_by(slug: params[:organism_slug])
    return head :not_found unless @organism

    @config = {
      "data_server": ENV['IGV_DATA_SERVER'],
      "organism_slug": @organism.slug,
      "organism_name": @organism.name
    }

    @track_data = @organism.lineage_variants_data(@config[:data_server], @config[:organism_slug])

    return unless @track_data[:data_fetched]

    @url_lineage = params[:lineage].presence
    @url_primer_sets = params[:primer_sets]&.split(',')

    @lineage_frequencies = @organism.recent_lineage_frequencies
  end
end
