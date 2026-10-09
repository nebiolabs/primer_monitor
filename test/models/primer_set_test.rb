# frozen_string_literal: true

require 'test_helper'

class PrimerSetTest < ActiveSupport::TestCase
  include ActionMailer::TestHelper

  test 'citation_url must be an http(s) link, since it is rendered as one' do
    primer_set = primer_sets(:one)

    primer_set.citation_url = 'javascript:alert(1)'

    assert_not primer_set.valid?
    assert_includes primer_set.errors[:citation_url], 'must start with http:// or https://'

    primer_set.citation_url = 'https://doi.org/10.1000/example'
    primer_set.valid?

    assert_empty primer_set.errors[:citation_url]
  end

  test 'admins are emailed only once the save has committed, so the mailer job can load the primer set' do
    primer_set = primer_sets(:one)

    PrimerSet.transaction do
      assert_no_enqueued_emails { primer_set.update!(name: 'renamed primer set') }
    end

    notified_admins = Role.find_by(name: 'administrator').users.where.not(email: ENV.fetch('ADMIN_EMAIL', nil))

    assert_enqueued_emails notified_admins.count
  end
end
