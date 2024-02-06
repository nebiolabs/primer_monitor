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

    @primer_sets = {}
    @default_tracks = []
    @lineage_sets = {}

    @data_fetched = false

    tracks_req = HTTParty.get("#{@config[:data_server]}/#{@config[:organism_slug]}/config/tracks.json")
    defaults_req = HTTParty.get("#{@config[:data_server]}/#{@config[:organism_slug]}/defaults.json")
    lineages_req = HTTParty.get("#{@config[:data_server]}/#{@config[:organism_slug]}/config/lineage_sets.json")

    return unless tracks_req.code == 200 && defaults_req.code == 200 && lineages_req.code == 200

    @data_fetched = true

    @primer_sets = JSON.parse(tracks_req.body)

    defaults_parsed = JSON.parse(defaults_req.body)
    @default_tracks = defaults_parsed['tracks']
    @default_lineage = defaults_parsed['lineage']

    if params.key? 'lineage'
      Rails.logger.warn params[:lineage]
      @default_lineage = params[:lineage]
    end

    if params.key? 'primer_sets'
      Rails.logger.warn params[:primer_sets]
      @default_tracks = params[:primer_sets].split(',')
    end

    @lineage_sets = JSON.parse(lineages_req.body)
  end
end
