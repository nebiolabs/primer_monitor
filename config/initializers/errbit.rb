Airbrake.configure do |config|
  config.host = 'http://langhorst-dev.neb.com'
  config.project_id = 1 # required, but any positive integer works
  config.project_key = 'bb47f2c39b1fbb060e913ccb5332805b'
  config.performance_stats = false

  # Uncomment for Rails apps
  config.environment = Rails.env
  config.ignore_environments = %w(development test)
end