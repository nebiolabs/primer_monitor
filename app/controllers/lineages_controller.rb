# frozen_string_literal: true

class LineagesController < ApplicationController
  require 'net/http'
  def index
    authorize! :index, LineagesController

    @config = {
      "data_server": ENV['IGV_DATA_SERVER'],
      "organism": 'SARS-CoV-2' # note: hardcoding SARS-CoV-2
    }

    tracks_url = URI("#{@config[:data_server]}/#{@config[:organism]}/misc/tracks.json")

    tracks = JSON.parse(Net::HTTP.get(tracks_url))

    @primer_sets = tracks['tracks']

    @default_tracks = tracks['default']

    lineages_url = URI("#{@config[:data_server]}/#{@config[:organism]}/misc/lineage_sets.json")

    @lineage_sets = JSON.parse(Net::HTTP.get(lineages_url))

    puts @default_tracks
  end
end
