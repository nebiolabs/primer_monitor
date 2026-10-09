# frozen_string_literal: true

require 'test_helper'

class ApplicationSystemTestCase < ActionDispatch::SystemTestCase
  driven_by :selenium, using: :headless_chrome, screen_size: [1400, 1400]

  # The igv.js pages load genomes and tracks from the IGV data server. Tests serve a small synthetic data set
  # (test/fixtures/files/igvstatic) from the app's own test server instead, so they need no network.
  IGV_FIXTURES = Rails.root.join('test/fixtures/files/igvstatic').to_s
  Capybara.app = Rack::Builder.new do
    map('/igvstatic') { run Rack::Files.new(IGV_FIXTURES) }
    run Rails.application
  end

  include Devise::Test::IntegrationHelpers

  setup do
    server = Capybara.current_session.server # Capybara picks a free port, so ask it which
    ENV['IGV_DATA_SERVER'] = "http://#{server.host}:#{server.port}/igvstatic"
  end
  teardown { ENV.delete('IGV_DATA_SERVER') }

  # igv.js 3 renders inside a shadow root on #igv, out of reach of normal selectors
  def igv_track_names
    page.evaluate_script(<<~JS)
      [...(document.getElementById('igv')?.shadowRoot?.querySelectorAll('.igv-track-label') || [])]
        .map(label => label.textContent.trim()).filter(Boolean)
    JS
  end

  # waits for igv to show exactly these tracks; order is ignored since primer set tracks load in parallel
  def assert_igv_tracks(expected, message = nil)
    nil
    deadline = Time.current + (Capybara.default_max_wait_time * 3)
    sleep 0.2 until (tracks = igv_track_names).sort == expected.sort || Time.current > deadline

    assert_equal expected.sort, tracks.sort,
                 message || "igv tracks: expected #{expected.inspect}, saw #{tracks.inspect}"
  end
end
