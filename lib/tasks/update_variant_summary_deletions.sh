
start=$1
end=$2
start=4101
end=21617

for blockstart in `seq $start 1000 $end`; do 
  blockend=$(echo "${blockstart} + 1000" | bc)
  echo psql -d primer_monitor -h ebase-db-c -c "\"with genpos as (select generate_series($blockstart,$blockend) pos)
update variant_summary_by_position vsp
set num_del = (select count(1) from variant_sites vs
where variant_type = 'D' and vsp.pos >= vs.ref_start and vsp.pos < vs.ref_end)
from genpos
where vsp.pos =genpos.pos;\""
done
