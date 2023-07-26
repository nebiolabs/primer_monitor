# frozen_string_literal: true

class LineageVariantsController < ApplicationController
  require 'net/http'
  def index
    authorize! :index, LineageVariantsController

    # NOTE: hardcoding SARS-CoV-2 taxid
    organism = Organism.find_by(ncbi_taxon_id: '2697049')

    @config = {
      "data_server": ENV['IGV_DATA_SERVER'],
      "organism_taxid": organism.ncbi_taxon_id,
      "organism_name": organism.name,
      "reference_accession": organism.reference_accession
    }

    tracks_url = URI("#{@config[:data_server]}/#{@config[:organism_taxid]}/config/tracks.json")

    @primer_sets = JSON.parse(Net::HTTP.get(tracks_url))

    defaults_url = URI("#{@config[:data_server]}/#{@config[:organism_taxid]}/defaults.json")

    defaults = JSON.parse(Net::HTTP.get(defaults_url))

    @default_tracks = defaults['tracks']

    lineages_url = URI("#{@config[:data_server]}/#{@config[:organism_taxid]}/config/lineage_sets.json")

    @lineage_sets = JSON.parse(Net::HTTP.get(lineages_url))

    puts @default_tracks
  end
end
