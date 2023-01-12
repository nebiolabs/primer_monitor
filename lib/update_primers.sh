#!/usr/bin/env bash

if (($# < 1)); then # if bowtie2 index missing
	echo "usage: get_primers.sh <bowtie2 index> [primer_set_ids...]" >&2
	exit 1
fi

if (($# < 2)); then # if primer set IDs missing
	echo "no primer set ids specified, exiting" >&2
	exit 1
fi

bt2_index=$1;
shift;

# you need to export DB_HOST, DB_NAME, and DB_USER before running this

db_csv=$(mktemp)


null_check="AND ref_start IS NULL";

if [ -n "$REALIGN" ]; then
  null_check="";
fi

for id in "$@"
do
    echo "Processing $id..." >&2;

    psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" -c "SELECT id, sequence FROM oligos WHERE primer_set_id=$id $null_check;" --csv -t | \
    awk 'BEGIN { FS="," }; {print ">" $1 "\n" $2}' | \
    bowtie2 -f --end-to-end --score-min L,-0.6,-1.5 -L 8 -x "$bt2_index" -U - | \
    samtools view -b | bedtools bamtobed -i - | awk '{print $4 "," $2 "," $3}' >> "$db_csv"

done

psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" >&2 <<CMDS
create temporary table tmp_oligo_positions (seq_id integer, ref_start integer, ref_end integer);
\copy tmp_oligo_positions from '$db_csv' with (format csv);
update oligos set ref_start=top.ref_start, ref_end=top.ref_end from tmp_oligo_positions top where oligos.id=top.seq_id;
CMDS

rm "$db_csv"