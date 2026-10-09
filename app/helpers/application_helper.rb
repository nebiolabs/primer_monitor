# frozen_string_literal: true

module ApplicationHelper
  # Links to a user-supplied URL only when it is http(s), so records saved before validation existed
  # cannot render javascript: or data: links. Anything else shows as plain text.
  def external_link(text, url)
    uri = URI.parse(url.to_s)
    return text unless uri.is_a?(URI::HTTP) && uri.host.present?

    link_to text, uri.to_s, rel: 'noopener', target: '_blank'
  rescue URI::InvalidURIError
    text
  end
end
