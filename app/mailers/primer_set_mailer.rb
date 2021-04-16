# frozen_string_literal: true

class PrimerSetMailer < ApplicationMailer
  default from: 'primer-monitor@neb.com'

  def updated_primer_set_email(email, primer_set)
    @url = primer_set_url(primer_set)
    mail(to: email, subject: 'Primer Set Created or Updated')
  end
end
