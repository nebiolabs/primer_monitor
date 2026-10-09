# frozen_string_literal: true

# Vendored stylesheets (vendor/assets/stylesheets): Bulma 1 and the DataTables / Tom Select CSS, which
# application.scss loads with @use (plain .css files are inlined when used without their extension).
# --quiet-deps hides deprecation warnings inside those vendored libraries (Bulma 1.0.4 still uses Sass's old if())
# while still reporting any in the app's own stylesheets.
Rails.application.config.dartsass.build_options |= %w[--load-path=vendor/assets/stylesheets --quiet-deps]
