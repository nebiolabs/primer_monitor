# frozen_string_literal: true

class PrimerSetsController < ApplicationController
  load_and_authorize_resource
  before_action :set_primer_set, only: %i[show edit update destroy]

  # GET /primer_sets
  # GET /primer_sets.json
  def index
    @primer_sets = PrimerSet.includes(:organism, :oligos).where(status: :complete)
    @subscriptions = PrimerSetSubscription.subscriptions_for_user_by_primer_set(current_user)
  end

  # GET /primer_sets/1
  # GET /primer_sets/1.json
  def show; end

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
        format.html do
          redirect_to @primer_set,
                      notice: 'Primer set was successfully added. It will be visible after processing and review.'
        end
        format.json { render :show, status: :created, location: @primer_set }
      else
        format.html { render :new }
        format.json { render json: @primer_set.errors, status: :unprocessable_entity }
      end
    end
  end

  # PATCH/PUT /primer_sets/1
  # PATCH/PUT /primer_sets/1.json
  def update
    respond_to do |format|
      if @primer_set.update(primer_set_params)
        format.html { redirect_to @primer_set, notice: 'Primer set was successfully updated.' }
        format.json { render :show, status: :ok, location: @primer_set }
      else
        format.html { render :edit }
        format.json { render json: @primer_set.errors, status: :unprocessable_entity }
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

  # Use callbacks to share common setup or constraints between actions.
  def set_primer_set
    @primer_set = PrimerSet.find(params[:id])
  end

  # Only allow a list of trusted parameters through.
  def primer_set_params
    params.require(:primer_set).permit(:name, :user_id, :amplification_method_id, :organism_id, :status,
                                       :citation_url, :doi,
                                       oligos_attributes: %i[id primer_set_id name short_name
                                                             locus category sequence _destroy])
  end
end
