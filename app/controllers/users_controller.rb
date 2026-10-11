# frozen_string_literal: true

class UsersController < ApplicationController
  load_and_authorize_resource

  before_action :set_user, only: %i[show edit update destroy]

  # GET /users
  # GET /users.json
  def index
    @users = User.all
  end

  # GET /users/1
  # GET /users/1.json
  def show; end

  # GET /users/1/edit
  def edit; end

  # PATCH/PUT /users/1
  # PATCH/PUT /users/1.json
  def update
    respond_to do |format|
      if @user.update(user_params)
        format.html { redirect_to @user, notice: 'User was successfully updated.' }
        format.json { render :show, status: :ok, location: @user }
      else
        format.html { render :edit, status: :unprocessable_content }
        format.json { render json: @user.errors, status: :unprocessable_content }
      end
    end
  end

  # DELETE /users/1
  # DELETE /users/1.json
  def destroy
    @user.destroy
    respond_to do |format|
      format.html { redirect_to users_url, notice: 'User was successfully destroyed.' }
      format.json { head :no_content }
    end
  end

  private

  def set_user
    @user = User.find(params.expect(:id))
  end

  # Only allow a list of trusted parameters through. Blank password fields keep the current password.
  # Only administrators can (de)activate a user.
  def user_params
    fields = %i[first last email password password_confirmation send_primer_updates]
    fields << :active if current_user.role?(:administrator)
    permitted = params.expect(user: fields)
    permitted[:password].blank? ? permitted.except(:password, :password_confirmation) : permitted
  end
end
