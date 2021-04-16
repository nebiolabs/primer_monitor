# frozen_string_literal: true

# Preview all emails at http://localhost:3000/rails/mailers/primer_set_mailer
class PrimerSetMailerPreview < ActionMailer::Preview
  def updated_primer_set_email
    PrimerSetMailer.updated_primer_set_email('me@example.com', PrimerSet.first)
  end

end
