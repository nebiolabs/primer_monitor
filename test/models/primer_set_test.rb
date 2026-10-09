# frozen_string_literal: true

require 'test_helper'

class PrimerSetTest < ActiveSupport::TestCase
  test 'citation_url must be an http(s) link, since it is rendered as one' do
    primer_set = primer_sets(:one)

    primer_set.citation_url = 'javascript:alert(1)'

    assert_not primer_set.valid?
    assert_includes primer_set.errors[:citation_url], 'must start with http:// or https://'

    primer_set.citation_url = 'https://doi.org/10.1000/example'
    primer_set.valid?

    assert_empty primer_set.errors[:citation_url]
  end
end
