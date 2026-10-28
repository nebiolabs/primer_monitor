# frozen_string_literal: true

# Geographic, lineage and primer summary of one reference position, linked from the notification emails. The page
# draws a map of the sequences and an igv.js diagram of the position; the .text format serves the diagram's tracks.
class VariantPositionsController < ApplicationController
  DEFAULT_DAYS = 30

  def show
    authorize! :show, VariantPositionsController

    @organism = Organism.find_by(slug: params[:organism_slug])
    return head :not_found unless @organism

    load_position

    respond_to do |format|
      format.html { @igv_config = igv_config }
      format.text do
        body = params[:track] == 'primers' ? @track.primer_bed : @track.variant_bed(params[:allele].to_i)
        body ? render(plain: body) : head(:not_found)
      end
    end
  end

  private

  def load_position
    @coordinate = params[:coordinate].to_i
    @days = params[:days].present? ? params[:days].to_i.clamp(1, 365) : DEFAULT_DAYS
    @summary = VariantPositions::Summary.call(organism: @organism, coordinate: @coordinate, days: @days)
    @track = VariantPositions::Track.new(organism: @organism, coordinate: @coordinate, days: @days)
  end

  def igv_config
    @organism.igv_config.merge(
      locus: @track.locus,
      variant_tracks: @track.alleles.each_with_index.map do |allele, index|
        { name: @track.allele_name(allele),
          url: organism_variant_position_url(@organism, coordinate: @coordinate, days: @days,
                                                        track: 'variant', allele: index, format: :text) }
      end,
      primer_track_url: organism_variant_position_url(@organism, coordinate: @coordinate, days: @days,
                                                                 track: 'primers', format: :text)
    )
  end
end
