# frozen_string_literal: true

class UserMailer < ApplicationMailer
  def verification_email
    @user = params[:user]
    mail(to: @user.email, subject: 'Primer Monitor Email Verification')
  end
end
