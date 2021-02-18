# frozen_string_literal: true

class ApplicationMailer < ActionMailer::Base
  default from: 'primer-monitor@neb.com'
  layout 'mailer'
end
