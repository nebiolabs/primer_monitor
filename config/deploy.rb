# config valid for current version and patch releases of Capistrano
lock '~> 3.15'

set :application, 'primer-monitor'
set :ssh_options, { forward_agent: true }
set :repo_url, 'git@github.com:bwlang/primer_monitor.git'
set :puma_service_name, 'puma'

set :whenever_roles, [:app, :backend]

# only migrate if there are migrations pending
set :conditionally_migrate, true

# Avoid recompiling all the assets on every deploy
# from https://coderwall.com/p/aridag/only-precompile-assets-when-necessary-rails-4-capistrano-3
# set the locations that we will look for changed assets to determine whether to precompile
set :assets_dependencies, %w[app/assets lib/assets vendor/assets Gemfile.lock config/routes.rb]

set :backend_deploy_to, ->{ fetch(:backend_deploy_path) }

# To get the backend path into whenever
set :whenever_environment, ->{ fetch(:rails_env, fetch(:stage, "production")) }
set :whenever_variables, ->{ "environment=#{fetch :whenever_environment}&backend_path=#{fetch(:backend_deploy_to)}" }
set :whenever_command, ->{ "bundle exec whenever -s #{fetch(:whenever_variables)}" }

# clear the previous precompile task
Rake::Task['deploy:assets:precompile'].clear_actions
class PrecompileRequired < StandardError
end

# Default branch is :master
ask :branch, `git rev-parse --abbrev-ref HEAD`.chomp

append :linked_files, 'config/database.yml', 'config/master.key', '.env'

append :linked_dirs,  'log', 'tmp/pids', 'tmp/cache', 'tmp/sockets', 'public/system'

set :keep_releases, 10

desc 'Deploy primer-monitor backend'
task :backend do
  invoke 'backend:git'
end

namespace :backend do
  desc 'Pull backend code from git'
  task :git do
    on roles(:backend) do
      within :backend_deploy_path do
        execute("cd #{fetch(:backend_deploy_path)} && git pull && git checkout #{fetch(:branch)}")
      end
    end
  end
end

namespace :deploy do
  desc 'Restart application'
  task :restart do
    on roles(:app), in: :sequence, wait: 5 do
      execute :touch, release_path.join('tmp/restart.txt')
    end
  end

  desc 'Run rake yarn:install'
  task :yarn_install do
    on roles(:app) do
      within release_path do
        execute("cd #{release_path} && yarn install && yarn upgrade")
      end
    end
  end

  desc 'Seed application'
  task :seed do
    on roles(:app) do
      with rails_env: fetch(:rails_env) do
        within release_path do
          execute :rake, 'db:seed'
        end
      end
    end
  end

  after :restart, :clear_cache do
    on roles(:app), in: :groups, limit: 3, wait: 10 do
      with rails_env: fetch(:rails_env) do
        within release_path do
          execute :rake, 'tmp:clear'
        end
      end
    end
  end

  after :restart, :restart_services do
    on roles(:app), in: :groups, limit: 3, wait: 10 do
      within release_path do
        execute 'mkdir -p tmp/sockets'
        execute "sudo /bin/systemctl restart #{fetch(:puma_service_name)}"
      end
    end
  end

  namespace :assets do
    desc 'Precompile assets'
    task :precompile do
      on roles(fetch(:assets_roles)) do
        within release_path do
          with rails_env: fetch(:rails_env) do
            # find the most recent release
            latest_release = capture(:ls, '-xr', releases_path).split[1]

            # precompile if this is the first deploy
            raise PrecompileRequired unless latest_release

            latest_release_path = releases_path.join(latest_release)

            # precompile if the previous deploy failed to finish precompiling
            begin
              execute(:ls, latest_release_path.join('assets_manifest_backup'))
            rescue StandardError
              raise(PrecompileRequired)
            end

            fetch(:assets_dependencies).each do |dep|
              # execute raises if there is a diff

              execute(:diff, '-Naur', release_path.join(dep),
                      latest_release_path.join(dep))
            rescue StandardError
              raise(PrecompileRequired)
            end

            info('Skipping asset precompile, no asset diff found')

            # copy over all of the assets from the last release
            execute(:cp, '-r', latest_release_path.join('public', fetch(:assets_prefix)),
                    release_path.join('public', fetch(:assets_prefix)))
          rescue PrecompileRequired
            execute(:rake, 'assets:precompile')
          end
        end
      end
    end
  end

  namespace :symlink do
    desc 'Link ref and defaults into igvstatic data'
    task :link_igvstatic do
      on roles :app do
        within release_path do
          execute <<~BASH
            for dir in #{release_path}/igvstatic/*; do \
              pushd $dir && \
              ln -sf #{shared_path}/igvstatic/$(basename \"$dir\")/* . && \
              popd
            done
          BASH
        end
      end
    end

    after :linked_files, :link_igvstatic

  end

  after :published, :restart
  #after 'deploy:restart_services', 'deploy:seed'
  after 'deploy:restart_services', 'backend'

end
