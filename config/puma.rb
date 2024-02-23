# frozen_string_literal: true

app_dir = File.expand_path('..', __dir__)
rails_env = ENV['RAILS_ENV'] || 'development'
environment rails_env

# Puma can serve each request in a thread from an internal thread pool.
# The `threads` method setting takes two numbers: a minimum and maximum.
# Any libraries that use thread pools should be configured to match
# the maximum value specified for Puma. Default is set to 5 threads for minimum
# and maximum; this matches the default thread size of Active Record.
#
max_threads_count = ENV.fetch('RAILS_MAX_THREADS', 5)
min_threads_count = ENV.fetch('RAILS_MIN_THREADS') { max_threads_count }
threads min_threads_count, max_threads_count

# Specifies the `environment` that Puma will run in.
if rails_env == 'development'
  workers 1
  # Specifies the `port` that Puma will listen on to receive requests; default is 3000.
  port ENV.fetch('PORT', 3000)
  # Specifies the `worker_timeout` threshold that Puma will use to wait before
  # terminating a worker in development environments.
  worker_timeout 3600
else
  workers 2
  bind "unix://#{app_dir}/tmp/sockets/puma.sock"
  stdout_redirect "#{app_dir}/log/puma.stdout.log", "#{app_dir}/log/puma.stderr.log", true
end
# Specifies the `pidfile` that Puma will use.
pidfile ENV.fetch('PIDFILE') { "#{app_dir}/tmp/pids/server.pid" }
state_path "#{app_dir}/tmp/pids/puma.state"

# Specifies the number of `workers` to boot in clustered mode.
# Workers are forked web server processes. If using threads and workers together
# the concurrency of the application would be max `threads` * `workers`.
# Workers do not work on JRuby or Windows (both of which do not support
# processes).
#
# workers ENV.fetch("WEB_CONCURRENCY") { 2 }

# Use the `preload_app!` method when specifying a `workers` number.
# This directive tells Puma to first boot the application and load code
# before forking the application. This takes advantage of Copy On Write
# process behavior so workers use less memory.
#
preload_app!

on_worker_boot do
  # Valid on Rails 4.1+ using the `config/database.yml` method of setting `pool` size
  ActiveRecord::Base.establish_connection
end

# Allow puma to be restarted by `rails restart` command.
plugin :tmp_restart
