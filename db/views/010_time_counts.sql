
CREATE MATERIALIZED VIEW public.time_counts AS
 WITH region_time_count AS (
         SELECT count(*) AS region_time_count,
            COALESCE(detailed_geo_locations.region, ''::character varying) AS region,
            COALESCE(fasta_records.date_collected, '1900-01-01'::date) AS date_collected
           FROM (public.fasta_records
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
          GROUP BY detailed_geo_locations.region, fasta_records.date_collected
        ), region_subregion_time_count AS (
         SELECT count(*) AS region_subregion_time_count,
            COALESCE(detailed_geo_locations.region, ''::character varying) AS region,
            COALESCE(detailed_geo_locations.subregion, ''::character varying) AS subregion,
            COALESCE(fasta_records.date_collected, '1900-01-01'::date) AS date_collected
           FROM (public.fasta_records
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, fasta_records.date_collected
        ), region_subregion_division_time_count AS (
         SELECT count(*) AS region_subregion_division_time_count,
            COALESCE(detailed_geo_locations.region, ''::character varying) AS region,
            COALESCE(detailed_geo_locations.subregion, ''::character varying) AS subregion,
            COALESCE(detailed_geo_locations.division, ''::character varying) AS division,
            COALESCE(fasta_records.date_collected, '1900-01-01'::date) AS date_collected
           FROM (public.fasta_records
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, fasta_records.date_collected
        ), region_subregion_division_subdivision_time_count AS (
         SELECT count(*) AS region_subregion_division_subdivision_time_count,
            COALESCE(detailed_geo_locations.region, ''::character varying) AS region,
            COALESCE(detailed_geo_locations.subregion, ''::character varying) AS subregion,
            COALESCE(detailed_geo_locations.division, ''::character varying) AS division,
            COALESCE(detailed_geo_locations.subdivision, ''::character varying) AS subdivision,
            COALESCE(fasta_records.date_collected, '1900-01-01'::date) AS date_collected
           FROM (public.fasta_records
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, detailed_geo_locations.subdivision, fasta_records.date_collected
        )
 SELECT region_subregion_division_subdivision_time_count.region,
    region_subregion_division_subdivision_time_count.subregion,
    region_subregion_division_subdivision_time_count.division,
    region_subregion_division_subdivision_time_count.subdivision,
    region_subregion_division_subdivision_time_count.date_collected,
    region_time_count.region_time_count,
    region_subregion_time_count.region_subregion_time_count,
    region_subregion_division_time_count.region_subregion_division_time_count,
    region_subregion_division_subdivision_time_count.region_subregion_division_subdivision_time_count
   FROM (((region_subregion_division_subdivision_time_count
     JOIN region_time_count ON ((((region_time_count.region)::text = (region_subregion_division_subdivision_time_count.region)::text) AND (region_time_count.date_collected = region_subregion_division_subdivision_time_count.date_collected))))
     JOIN region_subregion_time_count ON ((((region_subregion_time_count.region)::text = (region_subregion_division_subdivision_time_count.region)::text) AND ((region_subregion_time_count.subregion)::text = (region_subregion_division_subdivision_time_count.subregion)::text) AND (region_subregion_time_count.date_collected = region_subregion_division_subdivision_time_count.date_collected))))
     JOIN region_subregion_division_time_count ON ((((region_subregion_division_time_count.region)::text = (region_subregion_division_subdivision_time_count.region)::text) AND ((region_subregion_division_time_count.subregion)::text = (region_subregion_division_subdivision_time_count.subregion)::text) AND ((region_subregion_division_time_count.division)::text = (region_subregion_division_subdivision_time_count.division)::text) AND (region_subregion_division_time_count.date_collected = region_subregion_division_subdivision_time_count.date_collected))))
 ;

CREATE UNIQUE INDEX ON time_counts(region, subregion, division, subdivision, date_collected);

GRANT SELECT on time_counts to primer_monitor_ro;