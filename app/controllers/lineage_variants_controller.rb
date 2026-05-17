# frozen_string_literal: true

class LineageVariantsController < ApplicationController
  def index
    authorize! :index, LineageVariantsController

    @organism = Organism.find_by(slug: params[:organism_slug])
    return head :not_found unless @organism

    @config = {
      "data_server": ENV['IGV_DATA_SERVER'],
      "organism_slug": @organism.slug,
      "organism_name": @organism.name,
      "reference_accession": @organism.organism_taxa.first&.reference_accession
    }

    @track_data = @organism.lineage_variants_data(@config[:data_server], @config[:organism_slug])

    return unless @track_data[:data_fetched]

    @url_lineage = params[:lineage].presence
    @url_primer_sets = params[:primer_sets]&.split(',')

    prepare_lineage_display
    @config[:initial_lineage] = @effective_lineage
    @config[:initial_primer_sets] = @effective_primer_sets
  end

  def variant_overlaps
    authorize! :index, LineageVariantsController

    @organism = Organism.find_by(slug: params[:organism_slug])
    return head :not_found unless @organism

    lineage_param, primer_set_params = overlap_params
    return head :bad_request if invalid_overlap_params?(lineage_param, primer_set_params)

    variants = LineageVariants::OverlapResponseBuilder.call(
      organism: @organism,
      lineage_group_key: lineage_param,
      primer_set_names: primer_set_params
    )

    render json: { variants: }
  end

  private

  def overlap_params
    [params[:lineage].presence, Array(params[:primer_sets]).reject(&:blank?)]
  end

  def invalid_overlap_params?(lineage_param, primer_set_params)
    lineage_param.nil? || primer_set_params.empty?
  end

  def prepare_lineage_display
    state = LineageVariants::DisplayStateBuilder.call(
      organism: @organism,
      track_data: @track_data,
      url_lineage: @url_lineage,
      url_primer_sets: @url_primer_sets
    )

    @lineage_frequencies = state[:lineage_frequencies]
    @sorted_lineages = state[:sorted_lineages]
    @effective_lineage = state[:effective_lineage]
    @effective_primer_sets = state[:effective_primer_sets]
  end
end
