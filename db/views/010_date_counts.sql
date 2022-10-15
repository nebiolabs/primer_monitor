
CREATE VIEW public.date_counts AS
 SELECT count(*) AS total_count,
    fr.detailed_geo_location_id,
    fr.date_collected
   FROM public.fasta_records fr
   WHERE fr.date_collected > (select max(date_collected) - '24 weeks'::interval from fasta_records)
  GROUP BY fr.detailed_geo_location_id, fr.date_collected;

GRANT SELECT on date_counts to primer_monitor_ro;