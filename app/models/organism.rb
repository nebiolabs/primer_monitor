# frozen_string_literal: true

class Organism < ApplicationRecord
  has_many :blast_hits, dependent: :destroy
  has_many :primer_sets, dependent: :destroy
  has_many :organism_taxa, dependent: :destroy

  def to_s
    name
  end

  def to_param
    slug
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
    rescue JSON::ParserError, Errno::ECONNREFUSED, Errno::ENOENT, SocketError, Net::OpenTimeout, Net::ReadTimeout
      primer_sets = {}
    end

    [config, primer_sets]
  end

  def recent_lineage_frequencies(days: 30)
    most_recent = LineageInfo.where(organism_id: id).maximum(:last_seen)
    return {} unless most_recent

    LineageInfo
      .where(organism_id: id)
      .where('last_seen >= ?', most_recent - days)
      .pluck(:name, :times_seen)
      .to_h
  end

  def lineage_variants_data(data_server, organism_slug)
    tracks_req = HTTParty.get("#{data_server}/#{organism_slug}/config/tracks.json")
    defaults_req = HTTParty.get("#{data_server}/#{organism_slug}/defaults.json")
    lineages_req = HTTParty.get("#{data_server}/#{organism_slug}/config/lineage_sets.json")

    return { data_fetched: false } unless tracks_req.code == 200 && defaults_req.code == 200 && lineages_req.code == 200

    process_lineage_variants_data tracks_req.body, defaults_req.body, lineages_req.body
  rescue Errno::ECONNREFUSED, Errno::ENOENT, SocketError, Net::OpenTimeout, Net::ReadTimeout,
         HTTParty::Error, JSON::ParserError
    { data_fetched: false }
  end

  def process_lineage_variants_data(tracks_json, defaults_json, lineages_json)
    primer_sets = JSON.parse(tracks_json)

    defaults_parsed = JSON.parse(defaults_json)
    default_tracks = defaults_parsed['tracks']
    default_lineage = defaults_parsed['lineage']

    lineage_sets = JSON.parse(lineages_json)

    { data_fetched: true, primer_sets:, default_tracks:, lineage_sets:,
      default_lineage: }
  end
end
