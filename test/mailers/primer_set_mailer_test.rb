# frozen_string_literal: true

require 'test_helper'

class PrimerSetMailerTest < ActionMailer::TestCase
  test 'updated_primer_set_email' do
    # Create the email and store it for further assertions
    email = PrimerSetMailer.updated_primer_set_email('me@example.com', primer_sets(:one))

    # Send the email, then test that it got queued
    assert_emails 1 do
      email.deliver_now
    end

    # Test the body of the sent email contains what we expect it to
    assert_equal ['me@example.com'], email.to
    assert_equal ['primer-monitor@neb.com'], email.from
    assert_equal 'Primer Set Created or Updated', email.subject
  end
end
