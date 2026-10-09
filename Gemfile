# frozen_string_literal: true

source 'https://rubygems.org'
git_source(:github) { |repo| "https://github.com/#{repo}.git" }

ruby '3.4.2'
gem 'activerecord-import'

gem 'dotenv-rails'

gem 'httparty'

gem 'slop'

gem 'rails', '~> 8.1'

gem 'pg', '>=1.2.3'
# Use Puma as the app server
gem 'puma', '~> 8.0'
# Use SCSS for stylesheets
gem 'dartsass-rails'
# Build JSON APIs with ease. Read more: https://github.com/rails/jbuilder
gem 'jbuilder', '~> 2.7'
# Use Redis adapter to run Action Cable in production
# gem 'redis', '~> 4.0'
gem 'whenever'

# The original asset pipeline for Rails [https://github.com/rails/sprockets-rails]
gem 'sprockets-rails'

# Use JavaScript with ESM import maps [https://github.com/rails/importmap-rails]
gem 'importmap-rails'

# Hotwire's SPA-like page accelerator [https://turbo.hotwired.dev]
gem 'turbo-rails'

# Use Active Storage variant
# gem 'image_processing', '~> 1.2'

# Reduces boot times through caching; required in config/boot.rb
gem 'bootsnap', '>= 1.4.2', require: false

# for error reporting
gem 'airbrake'
# for authorization
gem 'cancancan', '~> 3.x'
# for authentication
gem 'devise', '~> 5.0'
gem 'omniauth-entra-id'
gem 'omniauth-google-oauth2'
gem 'omniauth-rails_csrf_protection', '~> 2.0'

# for nested form management
gem 'cocoon'

group :development, :test do
  # See https://guides.rubyonrails.org/debugging_rails_applications.html#debugging-with-the-debug-gem
  gem 'debug', platforms: %i[mri windows], require: 'debug/prelude'
  gem 'mocha'
end

group :development do
  gem 'bcrypt_pbkdf'
  gem 'brakeman', require: false
  gem 'bundler-audit', require: false
  gem 'capistrano', '~> 3.10', require: false
  gem 'capistrano-rails', '~> 1.3', require: false
  gem 'capistrano-rbenv', '~> 2.2'
  gem 'ed25519'
  gem 'listen', '~> 3.2'
  gem 'rubocop', '~> 1.91', require: false
  gem 'rubocop-capybara', '~> 3.0', require: false
  gem 'rubocop-minitest', '~> 0.41', require: false
  gem 'rubocop-performance', '~> 1.27', require: false
  gem 'rubocop-rails', '~> 2.38', require: false
  gem 'rubocop-rake', '~> 0.7', require: false

  # Access an interactive console on exception pages or by calling 'console' anywhere in the code.
  gem 'web-console', '>= 3.3.0'
  # for measuring test coverage
  gem 'simplecov'
  # for generating a search engine sitemap
  gem 'sitemap_generator'
  gem 'solargraph' # language server for code editors
  gem 'webrick'
end

group :test do
  # Adds support for Capybara system testing and selenium driver
  gem 'capybara', '>= 2.15'
  gem 'selenium-webdriver'
end

# Windows does not include zoneinfo files, so bundle the tzinfo-data gem
gem 'tzinfo-data', platforms: %i[windows jruby]
