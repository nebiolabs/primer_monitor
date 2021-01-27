# frozen_string_literal: true

class ApplicationMailer < ActionMailer::Base
  default from: 'primer-monitor-admin@neb.com'
  layout 'mailer'
end
