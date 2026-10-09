# frozen_string_literal: true

# Bulma 0.9 (vendor/assets/stylesheets/bulma) predates Sass modules, so silence its deprecation noise.
Rails.application.config.dartsass.build_options |= %w[
  --load-path=vendor/assets/stylesheets/bulma
  --quiet-deps
  --silence-deprecation=import
  --silence-deprecation=global-builtin
]
