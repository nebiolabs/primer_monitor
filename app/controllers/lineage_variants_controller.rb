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

    tracks_url = URI("#{@config[:data_server]}/#{@config[:organism_slug]}/config/tracks.json")

    @primer_sets = JSON.parse(Net::HTTP.get(tracks_url))

    defaults_url = URI("#{@config[:data_server]}/#{@config[:organism_slug]}/defaults.json")

    defaults = JSON.parse(Net::HTTP.get(defaults_url))

    @default_tracks = defaults['tracks']

    lineages_url = URI("#{@config[:data_server]}/#{@config[:organism_slug]}/config/lineage_sets.json")

    @lineage_sets = JSON.parse(Net::HTTP.get(lineages_url))
  end
end
