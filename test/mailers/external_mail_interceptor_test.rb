# frozen_string_literal: true

require 'test_helper'

class ExternalMailInterceptorTest < ActiveSupport::TestCase
  def build_mail(to:, cc: nil)
    mail = Mail.new(from: 'primer-monitor@neb.com', subject: 'Variants observed', to:)
    mail.cc = cc if cc
    mail
  end

  test 'keeps allowlisted recipients, including subdomains and display names' do
    mail = build_mail(to: ['Brad L <langhorst@neb.com>', 'someone@lab.NEB.com'])
    ExternalMailInterceptor.delivering_email(mail)

    assert_equal %w[langhorst@neb.com someone@lab.NEB.com], mail.to
    assert_includes mail[:to].to_s, 'Brad L <langhorst@neb.com>'
    assert mail.perform_deliveries
  end

  test 'drops external recipients and logs each one' do
    mail = build_mail(to: ['langhorst@neb.com', 'anastacia.griego@state.nm.us'], cc: 'x@notneb.com')
    logged = capture_warnings { ExternalMailInterceptor.delivering_email(mail) }

    assert_equal ['langhorst@neb.com'], mail.to
    assert_nil mail.cc
    assert mail.perform_deliveries
    assert_match(/Dropped external mail to anastacia\.griego@state\.nm\.us \(to\): Variants observed/, logged)
    assert_match(/Dropped external mail to x@notneb\.com \(cc\)/, logged)
  end

  test 'cancels a message with only external recipients' do
    mail = build_mail(to: 'anastacia.griego@state.nm.us')
    ExternalMailInterceptor.delivering_email(mail)

    assert_not mail.perform_deliveries
  end

  test 'MAIL_ALLOWED_DOMAINS replaces the default allowlist' do
    ENV['MAIL_ALLOWED_DOMAINS'] = 'example.org, state.nm.us'

    assert ExternalMailInterceptor.allowed?('a@state.nm.us')
    assert ExternalMailInterceptor.allowed?('b@example.org')
    assert_not ExternalMailInterceptor.allowed?('c@neb.com')
  ensure
    ENV.delete('MAIL_ALLOWED_DOMAINS')
  end

  private

  def capture_warnings
    output = StringIO.new
    original = Rails.logger
    Rails.logger = ActiveSupport::Logger.new(output)
    yield
    output.string
  ensure
    Rails.logger = original
  end
end
