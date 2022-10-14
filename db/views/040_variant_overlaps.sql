CREATE MATERIALIZED VIEW public.variant_overlaps AS
SELECT oligos.id                                                           AS oligo_id,
       oligos.name                                                         AS oligo_name,
       oligos.ref_start                                                    AS oligo_start,
       oligos.ref_end                                                      AS oligo_end,
       oligos.short_name                                                   AS oligo_short_name,
       oligos.locus                                                        AS oligo_locus_name,
       oligos.category                                                     AS oligo_primer_type,
       primer_sets.id                                                      AS primer_set_id,
       primer_sets.name                                                    AS primer_set_name,
       variant_sites.id                                                    AS variant_id,
       variant_sites.variant_type,
       variant_sites.variant,
       variant_sites.ref_start                                             AS variant_start,
       variant_sites.ref_end                                               AS variant_end,
       COALESCE(detailed_geo_locations.region, ''::character varying)      AS region,
       COALESCE(detailed_geo_locations.subregion, ''::character varying)   AS subregion,
       COALESCE(detailed_geo_locations.division, ''::character varying)    AS division,
       COALESCE(detailed_geo_locations.subdivision, ''::character varying) AS subdivision,
       detailed_geo_locations.id                                           AS detailed_geo_location_id,
       COALESCE(fasta_records.date_collected, '1900-01-01'::date)          AS date_collected
FROM ((((public.variant_sites
    JOIN public.fasta_records ON (variant_sites.fasta_record_id = fasta_records.id
        AND date_collected > (select max(date_collected) - '24 weeks'::interval from fasta_records fr)
        )
    JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
    JOIN public.oligos ON ((NOT ((oligos.ref_start >= variant_sites.ref_end) OR
                                 (oligos.ref_end <= variant_sites.ref_start)))))
    JOIN public.primer_sets ON (oligos.primer_set_id = primer_sets.id AND primer_sets.status = 'complete')))
WHERE (variant_sites.usable_del_or_snp = true)
UNION ALL
SELECT oligos.id                                                           AS oligo_id,
       oligos.name                                                         AS oligo_name,
       oligos.ref_start                                                    AS oligo_start,
       oligos.ref_end                                                      AS oligo_end,
       oligos.short_name                                                   AS oligo_short_name,
       oligos.locus                                                        AS oligo_locus_name,
       oligos.category                                                     AS oligo_primer_type,
       primer_sets.id                                                      AS primer_set_id,
       primer_sets.name                                                    AS primer_set_name,
       variant_sites.id                                                    AS variant_id,
       variant_sites.variant_type,
       variant_sites.variant,
       variant_sites.ref_start                                             AS variant_start,
       variant_sites.ref_end                                               AS variant_end,
       COALESCE(detailed_geo_locations.region, ''::character varying)      AS region,
       COALESCE(detailed_geo_locations.subregion, ''::character varying)   AS subregion,
       COALESCE(detailed_geo_locations.division, ''::character varying)    AS division,
       COALESCE(detailed_geo_locations.subdivision, ''::character varying) AS subdivision,
       detailed_geo_locations.id                                           AS detailed_geo_location_id,
       COALESCE(fasta_records.date_collected, '1900-01-01'::date)          AS date_collected
FROM ((((public.variant_sites
    JOIN public.fasta_records ON (
                variant_sites.fasta_record_id = fasta_records.id
            AND date_collected >
                (select max(date_collected) - '24 weeks'::interval
                 from fasta_records fr))
    JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
    JOIN public.oligos ON (((variant_sites.ref_start > oligos.ref_start) AND
                            (variant_sites.ref_start <= oligos.ref_end))))
    JOIN public.primer_sets ON ((oligos.primer_set_id = primer_sets.id)) AND primer_sets.status = 'complete'))
WHERE (variant_sites.usable_insertion = true)
WITH NO DATA
;
CREATE index on variant_overlaps (variant_type);
GRANT SELECT on variant_overlaps to primer_monitor_ro;
ALTER MATERIALIZED VIEW variant_overlaps owner to primer_monitor;
