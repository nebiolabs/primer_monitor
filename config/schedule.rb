# Use this file to easily define all of your cron jobs.
#
# It's helpful, but not entirely necessary to understand cron before proceeding.
# http://en.wikipedia.org/wiki/Cron

# Example:
#
# set :output, "/path/to/my/cron_log.log"
#
# every 2.hours do
#   command "/usr/bin/some_great_command"
#   runner "MyModel.some_method"
#   rake "some:great:rake:task"
# end
#
# every 4.days do
#   runner "AnotherModel.prune_old_records"
# end

# Learn more: http://github.com/javan/whenever

set :output, "/var/www/primer-monitor/shared/log/cron.log"

every 1.day, at: ['7:00 am'], roles: [:app] do
  rake 'notifications:send'
end

every 1.day, at: ['2:00 am'], roles: [:backend] do
  command "source #{backend_path}/.env && export SGE_QMASTER_PORT && export SGE_ROOT && \
  /usr/bin/qsub -S /bin/bash -v \"NXF_JAVA_HOME=$NXF_JAVA_HOME\" \
  #{backend_path}/lib/backend_scripts/primer-monitor.sh", output: "\"$BACKEND_SCRATCH_PATH/download_cron.log\""
end

every 2.weeks, at: ['3:00 am'], roles: [:backend] do
  command "source #{backend_path}/.env && export SGE_QMASTER_PORT && export SGE_ROOT && \
  /usr/bin/qsub -S /bin/bash -v \"NXF_JAVA_HOME=$NXF_JAVA_HOME\" \
  #{backend_path}/lib/backend_scripts/recall-pangolin.sh", output: "\"$BACKEND_SCRATCH_PATH/pangolin_cron.log\""
end

every 1.minute, roles: [:backend] do
  command "source #{backend_path}/.env && export SGE_QMASTER_PORT && export SGE_ROOT && \
  /usr/bin/qsub -S /bin/bash -v \"NXF_JAVA_HOME=$NXF_JAVA_HOME\" \
  #{backend_path}/lib/backend_scripts/process-primer-sets.sh", output: "\"$BACKEND_SCRATCH_PATH/new_primers_cron.log\""
end
