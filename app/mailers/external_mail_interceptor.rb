# frozen_string_literal: true

# Outside production the database holds copies of real users, so staging and development must never email them.
# Removes every recipient outside MAIL_ALLOWED_DOMAINS (comma-separated, default neb.com; subdomains count),
# logging each one, and cancels the message if no recipient is left.
class ExternalMailInterceptor
  FIELDS = %i[to cc bcc].freeze

  def self.active?
    Array(Rails.application.config.action_mailer.interceptors).map(&:to_s).include?(name)
  end

  def self.allowed_domains
    ENV.fetch('MAIL_ALLOWED_DOMAINS', 'neb.com').split(',').map { |domain| domain.strip.downcase }.compact_blank
  end

  def self.allowed?(address)
    domain = address.to_s.split('@').last.to_s.downcase
    allowed_domains.any? { |allowed| domain == allowed || domain.end_with?(".#{allowed}") }
  end

  def self.delivering_email(message)
    FIELDS.each do |field|
      next if message[field].nil?

      recipients = message[field].addrs
      dropped, kept = recipients.partition { |addr| !allowed?(addr.address) }
      dropped.each do |addr|
        Rails.logger.warn("Dropped external mail to #{addr.address} (#{field}): #{message.subject}")
      end
      message[field] = kept.empty? ? nil : kept.map(&:format) if dropped.any?
    end
    message.perform_deliveries = false if FIELDS.all? { |field| message[field].nil? }
  end
end
