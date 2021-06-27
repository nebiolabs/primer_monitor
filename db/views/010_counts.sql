
CREATE MATERIALIZED VIEW public.counts AS
 WITH region_count AS (
         SELECT count(*) AS region_count,
            COALESCE(detailed_geo_locations.region, ''::character varying) AS region
           FROM (public.fasta_records
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
          GROUP BY detailed_geo_locations.region
        ), region_subregion_count AS (
         SELECT count(*) AS region_subregion_count,
            COALESCE(detailed_geo_locations.region, ''::character varying) AS region,
            COALESCE(detailed_geo_locations.subregion, ''::character varying) AS subregion
           FROM (public.fasta_records
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion
        ), region_subregion_division_count AS (
         SELECT count(*) AS region_subregion_division_count,
            COALESCE(detailed_geo_locations.region, ''::character varying) AS region,
            COALESCE(detailed_geo_locations.subregion, ''::character varying) AS subregion,
            COALESCE(detailed_geo_locations.division, ''::character varying) AS division
           FROM (public.fasta_records
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division
        ), region_subregion_division_subdivision_count AS (
         SELECT count(*) AS region_subregion_division_subdivision_count,
            COALESCE(detailed_geo_locations.region, ''::character varying) AS region,
            COALESCE(detailed_geo_locations.subregion, ''::character varying) AS subregion,
            COALESCE(detailed_geo_locations.division, ''::character varying) AS division,
            COALESCE(detailed_geo_locations.subdivision, ''::character varying) AS subdivision
           FROM (public.fasta_records
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, detailed_geo_locations.subdivision
        )
 SELECT region_subregion_division_subdivision_count.region,
    region_subregion_division_subdivision_count.subregion,
    region_subregion_division_subdivision_count.division,
    region_subregion_division_subdivision_count.subdivision,
    region_count.region_count,
    region_subregion_count.region_subregion_count,
    region_subregion_division_count.region_subregion_division_count,
    region_subregion_division_subdivision_count.region_subregion_division_subdivision_count
   FROM (((region_subregion_division_subdivision_count
     JOIN region_count ON (((region_subregion_division_subdivision_count.region)::text = (region_count.region)::text)))
     JOIN region_subregion_count ON ((((region_subregion_count.region)::text = (region_subregion_division_subdivision_count.region)::text) AND ((region_subregion_count.subregion)::text = (region_subregion_division_subdivision_count.subregion)::text))))
     JOIN region_subregion_division_count ON ((((region_subregion_division_count.region)::text = (region_subregion_division_subdivision_count.region)::text) AND ((region_subregion_division_count.subregion)::text = (region_subregion_division_subdivision_count.subregion)::text) AND ((region_subregion_division_count.division)::text = (region_subregion_division_subdivision_count.division)::text))))
;
CREATE UNIQUE INDEX ON counts(region, subregion, division, subdivision);
