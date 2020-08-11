class PrimerSetsController < ApplicationController
  before_action :set_primer_set, only: [:show, :edit, :update, :destroy]

  # GET /primer_sets
  # GET /primer_sets.json
  def index
    @primer_sets = PrimerSet.all
  end

  # GET /primer_sets/1
  # GET /primer_sets/1.json
  def show
  end

  # GET /primer_sets/new
  def new
    @primer_set = PrimerSet.new
  end

  # GET /primer_sets/1/edit
  def edit
  end

  # POST /primer_sets
  # POST /primer_sets.json
  def create
    @primer_set = PrimerSet.new(primer_set_params)

    respond_to do |format|
      if @primer_set.save
        format.html { redirect_to @primer_set, notice: 'Primer set was successfully created.' }
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
      params.require(:primer_set).permit(:name, :creator_id)
    end
end
