# frozen_string_literal: true

class LineagesController < ApplicationController
  require 'net/http'
  def index
    authorize! :index, LineagesController

    @config = { "data_server": ENV['IGV_DATA_SERVER'] }

    tracks_url = URI("#{@config[:data_server]}/misc/tracks.json")

    @primer_sets = JSON.parse(Net::HTTP.get(tracks_url))

    lineages_url = URI("#{@config[:data_server]}/misc/lineage_sets.json")

    @lineage_sets = JSON.parse(Net::HTTP.get(lineages_url))
  end
end
