# frozen_string_literal: true

require 'test_helper'

class ApplicationHelperTest < ActionView::TestCase
  test 'external_link links http(s) URLs' do
    assert_dom_equal '<a rel="noopener" target="_blank" href="https://doi.org/10.1/x">paper</a>',
                     external_link('paper', 'https://doi.org/10.1/x')
  end

  test 'external_link shows anything else as text' do
    ['javascript:alert(1)', 'data:text/html,hi', 'not a url', nil].each do |url|
      assert_equal 'paper', external_link('paper', url), url.inspect
    end
  end
end
