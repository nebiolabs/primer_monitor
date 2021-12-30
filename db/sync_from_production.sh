#!/bin/sh
if [ ! -e /var/run/postgresql ]; then
    export PGHOST=/tmp
fi

default_db=primer_monitor_dev
default_user=`whoami`

DBNAME=${1:-$default_db}
USER=${2:-$default_user}

echo "replacing $DBNAME with with production snapshot"
rsync -P ebase-db-c:/var/backups/postgres/primer_monitor.pgdump /tmp


echo "SELECT
    pg_terminate_backend(pid)
FROM
    pg_stat_activity
WHERE -- don't kill my own connection!
    datname = '${DBNAME}' and
    pid <> pg_backend_pid(); " | psql -U $USER $DBNAME

echo "drop db"
dropdb -U $USER $DBNAME  || (echo 'failed to drop db' && exit 1)

echo "create db"
createdb -U $USER $DBNAME
pg_restore -j 8 -U $USER -d $DBNAME /tmp/primer_monitor.pgdump

bundle exec rake db:seed
bundle exec rake db:migrate
bundle exec rake db:seed
