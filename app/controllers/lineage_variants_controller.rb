# frozen_string_literal: true

class LineageVariantsController < ApplicationController
  require 'net/http'
  def index
    authorize! :index, LineageVariantsController

    organism = Organism.find_by(slug: params[:organism_name])

    @config = {
      "data_server": ENV['IGV_DATA_SERVER'],
      "organism_slug": organism.slug,
      "organism_name": organism.name
    }

    @primer_sets = {}
    @default_tracks = []
    @lineage_sets = {}

    tracks_url = URI("#{@config[:data_server]}/#{@config[:organism_slug]}/config/tracks.json")
    tracks_data = Net::HTTP.get_response(tracks_url)
    @primer_sets = JSON.parse(tracks_data.body) if tracks_data.code == '200'

    defaults_url = URI("#{@config[:data_server]}/#{@config[:organism_slug]}/defaults.json")
    defaults_data = Net::HTTP.get_response(defaults_url)
    defaults_parsed = JSON.parse(defaults_data.body)
    @default_tracks = defaults_parsed['tracks'] if defaults_data.code == '200'
    @default_lineage = defaults_parsed['lineage'] if defaults_data.code == '200'

    if params.key? 'lineage'
      Rails.logger.warn params[:lineage]
      @default_lineage = params[:lineage]
    end

    if params.key? 'primer_sets'
      Rails.logger.warn params[:primer_sets]
      @default_tracks = params[:primer_sets].split(',')
    end


    lineages_url = URI("#{@config[:data_server]}/#{@config[:organism_slug]}/config/lineage_sets.json")
    lineages_data = Net::HTTP.get_response(lineages_url)
    @lineage_sets = JSON.parse(lineages_data.body) if lineages_data.code == '200'

  end
end
