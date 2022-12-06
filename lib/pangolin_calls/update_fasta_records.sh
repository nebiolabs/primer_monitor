if (($# < 1)); then
    echo "usage: ./update_fasta_records.sh <pangolin CSV path>" >&2;
    exit 1;
fi

psql -h $DB_HOST -d $DB_NAME -U $DB_USER -f update_fasta_records.sql -v pangolin_csv_path=\'$1\';