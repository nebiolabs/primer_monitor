# frozen_string_literal: true

class LineagesController < ApplicationController
  DEFAULT_DAYS = 90
  MAX_DAYS = 365

  def index
    authorize! :index, LineagesController
    if params.key? :organism_slug
      @organism = Organism.find_by(slug: params[:organism_slug])

      @lineages = LineageInfo.where(organism_id: @organism.id).select(:name, :times_seen, :last_seen).to_a
    else
      # hardcoded legacy redirect
      redirect_to organism_lineages_url(Organism.find_by(slug: 'sars-cov-2').slug), status: :moved_permanently
    end
  end

  def show
    authorize! :show, LineagesController
    @organism = Organism.find_by(slug: params[:organism_slug])
    @lineage = Lineage.find_by(name: params[:name])
    @lineage_info = LineageInfo.find_by(organism_id: @organism.id, name: @lineage.name)

    load_window
    @geography = Lineages::Geography.call(organism: @organism, lineage: @lineage, since: @since)
    @impacts = Lineages::PrimerImpacts.new(organism: @organism, lineage: @lineage, since: @since)
    @igv_config = @organism.igv_config.merge(primer_set_tracks: primer_set_tracks)
  end

  # one primer set's primers as a BED track, for the lineage page's diagram
  def primer_track
    authorize! :show, LineagesController
    @organism = Organism.find_by(slug: params[:organism_slug])
    @lineage = Lineage.find_by(name: params[:lineage], organism_id: @organism&.id)
    return head :not_found if @lineage.nil?

    load_window
    impacts = Lineages::PrimerImpacts.new(organism: @organism, lineage: @lineage, since: @since)
    primer_set = impacts.set(params[:set].to_i)
    return head :not_found if primer_set.nil?

    render plain: impacts.bed(primer_set)
  end

  def to_param
    name
  end

  private

  def load_window
    @days = params[:days].present? ? params[:days].to_i.clamp(1, MAX_DAYS) : DEFAULT_DAYS
    @since = Date.current - @days
  end

  def primer_set_tracks
    @impacts.sets.map do |primer_set|
      { name: "#{primer_set.name}, impact #{primer_set.impact ? format('%.1f', primer_set.impact) : 'none'}",
        color: VariantPositions::Disruption.color(primer_set.impact),
        url: organism_lineage_primer_track_url(@organism, lineage: @lineage.name, set: primer_set.id, days: @days,
                                                          format: :text) }
    end
  end
end
