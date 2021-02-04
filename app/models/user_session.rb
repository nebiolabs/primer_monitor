# frozen_string_literal: true

class UserSession < Authlogic::Session::Base

  validate :check_if_verified



  private

  def check_if_verified
    errors.add(:base, "Account is not yet verified") unless attempted_record && attempted_record.confirmed
  end

end
