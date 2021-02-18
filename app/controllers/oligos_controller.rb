# frozen_string_literal: true

class OligosController < ApplicationController
  load_and_authorize_resource
  before_action :set_oligo, only: %i[show edit update destroy]

  # GET /oligos
  # GET /oligos.json
  def index
    @oligos = Oligo.all
  end

  # GET /oligos/1
  # GET /oligos/1.json
  def show; end

  # GET /oligos/new
  def new
    @oligo = Oligo.new
  end

  # GET /oligos/1/edit
  def edit; end

  # POST /oligos
  # POST /oligos.json
  def create
    @oligo = Oligo.new(oligo_params)

    respond_to do |format|
      if @oligo.save
        format.html { redirect_to @oligo, notice: 'Oligo was successfully created.' }
        format.json { render :show, status: :created, location: @oligo }
      else
        format.html { render :new }
        format.json { render json: @oligo.errors, status: :unprocessable_entity }
      end
    end
  end

  # PATCH/PUT /oligos/1
  # PATCH/PUT /oligos/1.json
  def update
    respond_to do |format|
      if @oligo.update(oligo_params)
        format.html { redirect_to @oligo, notice: 'Oligo was successfully updated.' }
        format.json { render :show, status: :ok, location: @oligo }
      else
        format.html { render :edit }
        format.json { render json: @oligo.errors, status: :unprocessable_entity }
      end
    end
  end

  # DELETE /oligos/1
  # DELETE /oligos/1.json
  def destroy
    @oligo.destroy
    respond_to do |format|
      format.html { redirect_to oligos_url, notice: 'Oligo was successfully destroyed.' }
      format.json { head :no_content }
    end
  end

  private

  # Use callbacks to share common setup or constraints between actions.
  def set_oligo
    @oligo = Oligo.find(params[:id])
  end

  # Only allow a list of trusted parameters through.
  def oligo_params
    params.require(:oligo).permit(:name, :sequence)
  end
end
