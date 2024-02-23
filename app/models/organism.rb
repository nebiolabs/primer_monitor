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
end
