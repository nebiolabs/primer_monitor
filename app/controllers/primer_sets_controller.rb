# frozen_string_literal: true

class PrimerSetsController < ApplicationController
  load_and_authorize_resource
  before_action :set_primer_set, only: %i[show edit update destroy]

  # GET /primer_sets
  # GET /primer_sets.json
  def index
    @organism = Organism.find_by(slug: params[:organism_slug])

    @subscriptions = PrimerSetSubscription.subscriptions_for_user_by_primer_set(current_user)

    @config, @primer_sets = @organism.primer_sets_config
  end

  # GET /primer_sets/1
  # GET /primer_sets/1.json
  def show
    @config, primer_sets = @primer_set.organism.primer_sets_config

    return unless primer_sets.key? @primer_set.name

    @config['primer_set_display_name'] = @primer_set.name
    @config['primer_set_name'] = primer_sets[@primer_set.name]
    prepare_lineage_display
  end

  # GET /primer_sets/new
  def new
    @primer_set = PrimerSet.new
  end

  # GET /primer_sets/1/edit
  def edit; end

  # POST /primer_sets
  # POST /primer_sets.json
  def create
    @primer_set = PrimerSet.new(primer_set_params)
    @primer_set.user ||= current_user

    respond_to do |format|
      if @primer_set.save
        # start align script in the background
        @primer_set.align_primers
        format.html do
          redirect_to edit_primer_set_url(@primer_set),
                      notice: 'Primer set was successfully added. It will be visible after processing and review.'
        end
        format.json { render :show, status: :created, location: @primer_set }
      else
        format.html { render :new, status: :unprocessable_content }
        format.json { render json: @primer_set.errors, status: :unprocessable_content }
      end
    end
  end

  # PATCH/PUT /primer_sets/1
  # PATCH/PUT /primer_sets/1.json
  def update
    respond_to do |format|
      if @primer_set.update(primer_set_params)
        @primer_set.align_primers
        format.html { redirect_to edit_primer_set_url(@primer_set), notice: 'Primer set was successfully updated.' }
        format.json { render :show, status: :ok, location: @primer_set }
      else
        format.html { render :edit, status: :unprocessable_content }
        format.json { render json: @primer_set.errors, status: :unprocessable_content }
      end
    end
  end

  # DELETE /primer_sets/1
  # DELETE /primer_sets/1.json
  def destroy
    @primer_set.destroy
    respond_to do |format|
      format.html { redirect_to primer_sets_url, notice: 'Primer set was successfully destroyed.' }
      format.json { head :no_content }
    end
  end

  private

  # picks the lineage whose variants are shown against this primer set: the one in the URL, else the most
  # prevalent recent one. Without the lineage data the page falls back to showing only the primer positions.
  def prepare_lineage_display
    organism = @primer_set.organism
    track_data = organism.lineage_variants_data(@config[:data_server], @config[:organism_slug])
    return unless track_data[:data_fetched]

    state = LineageVariants::DisplayStateBuilder.call(
      organism:, track_data:, url_lineage: params[:lineage].presence, url_primer_sets: [@config['primer_set_name']]
    )
    @lineage_sets = track_data[:lineage_sets]
    @lineage_frequencies = state[:lineage_frequencies]
    @sorted_lineages = state[:sorted_lineages]
    @config[:initial_lineage] = state[:effective_lineage]
    @config[:lineage_names] = @lineage_sets
    @config[:variant_overlaps_path] = variant_overlaps_organism_lineage_variants_path(organism, format: :json)
  end

  # Use callbacks to share common setup or constraints between actions.
  def set_primer_set
    @primer_set = PrimerSet.find(params.expect(:id))
  end

  # Only allow a list of trusted parameters through.
  def primer_set_params
    # [[...]] accepts the index-keyed hash that fields_for submits for nested oligos
    params.expect(primer_set: [:name, :user_id, :amplification_method_id, :organism_id, :status,
                               :citation_url, :doi,
                               { oligos_attributes: [%i[id primer_set_id name short_name
                                                        locus category sequence _destroy]] }])
  end
end
