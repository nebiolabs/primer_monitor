Rake::Task["db:setup"].enhance do
  Rake::Task["db:ensure_db_roles_exist"].invoke
end

Rake::Task["db:reset"].enhance do
  Rake::Task["db:ensure_db_roles_exist"].invoke
end

namespace :db do
  desc "Create database roles"
  task ensure_db_roles_exist: :environment do
    %w[primer_monitor primer_monitor_ro].each do |r|
      role_exists = ActiveRecord::Base.connection.execute("SELECT 1 FROM pg_roles WHERE rolname = '#{r}'").any?
      if role_exists
        puts "Role #{r} already exists"
      else
        puts "Creating role #{r}"
        ActiveRecord::Base.connection.execute <<-SQL
          CREATE ROLE #{r};
        SQL
      end
    end
  end

  desc "Drop database roles"
  task drop_roles: :environment do
    ActiveRecord::Base.connection.execute <<-SQL
      DROP ROLE IF EXISTS primer_monitor;
      DROP ROLE IF EXISTS primer_monitor_ro;
    SQL
  end
end
