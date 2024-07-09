# frozen_string_literal: true

class Organism < ApplicationRecord
  has_many :blast_hits, dependent: :destroy
  has_many :primer_sets, dependent: :destroy
  has_many :organism_taxa, dependent: :destroy

  def to_s
    name
  end

  def full_name
    name + (self.alias.blank? ? '' : " (#{self.alias})")
  end

  def primer_sets_config
    config = {
      "data_server": ENV['IGV_DATA_SERVER'],
      "organism_slug": slug,
      "organism_name": name
    }

    tracks_url = URI("#{config[:data_server]}/#{config[:organism_slug]}/config/tracks.json")

    begin
      primer_sets = JSON.parse(Net::HTTP.get(tracks_url)).invert
    rescue JSON::ParserError
      # cannot parse JSON properly, treat it as empty
      primer_sets = {}
    end

    [config, primer_sets]
  end

  def lineage_variants_data(data_server, organism_slug)
    tracks_req = HTTParty.get("#{data_server}/#{organism_slug}/config/tracks.json")
    defaults_req = HTTParty.get("#{data_server}/#{organism_slug}/defaults.json")
    lineages_req = HTTParty.get("#{data_server}/#{organism_slug}/config/lineage_sets.json")

    unless tracks_req.code == 200 && defaults_req.code == 200 && lineages_req.code == 200
      # data_fetched, @primer_sets, @default_tracks, @lineage_sets, @default_lineage
      return { data_fetched: false }
    end

    primer_sets = JSON.parse(tracks_req.body)

    defaults_parsed = JSON.parse(defaults_req.body)
    default_tracks = defaults_parsed['tracks']
    default_lineage = defaults_parsed['lineage']

    lineage_sets = JSON.parse(lineages_req.body)

    { data_fetched: true, primer_sets: primer_sets, default_tracks: default_tracks, lineage_sets: lineage_sets,
      default_lineage: default_lineage }
  end
end
