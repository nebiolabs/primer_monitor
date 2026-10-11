CREATE MATERIALIZED VIEW public.oligo_variant_overlaps AS
SELECT variant_overlaps.oligo_id,
       variant_overlaps.oligo_name,
       variant_overlaps.oligo_start,
       variant_overlaps.oligo_end,
       variant_overlaps.oligo_short_name,
       variant_overlaps.oligo_locus_name,
       variant_overlaps.oligo_primer_type,
       variant_overlaps.primer_set_id,
       variant_overlaps.primer_set_name,
       variant_overlaps.variant_id,
       variant_overlaps.variant_type,
       variant_overlaps.variant,
       variant_overlaps.variant_start,
       variant_overlaps.variant_end,
       variant_overlaps.region,
       variant_overlaps.subregion,
       variant_overlaps.division,
       variant_overlaps.subdivision,
       variant_overlaps.detailed_geo_location_id,
       variant_overlaps.date_collected,
       generate_series(lower((numrange((variant_overlaps.oligo_start)::numeric, (variant_overlaps.oligo_end)::numeric) *
                              numrange((variant_overlaps.variant_start)::numeric,
                                       (variant_overlaps.variant_end)::numeric))),
                       (upper((numrange((variant_overlaps.oligo_start)::numeric,
                                        (variant_overlaps.oligo_end)::numeric) *
                               numrange((variant_overlaps.variant_start)::numeric,
                                        (variant_overlaps.variant_end)::numeric))) - (1)::numeric)) AS coords
FROM public.variant_overlaps variant_overlaps
WHERE ((variant_overlaps.variant_type)::text <> 'I'::text)
UNION ALL
SELECT variant_overlaps.oligo_id,
       variant_overlaps.oligo_name,
       variant_overlaps.oligo_start,
       variant_overlaps.oligo_end,
       variant_overlaps.oligo_short_name,
       variant_overlaps.oligo_locus_name,
       variant_overlaps.oligo_primer_type,
       variant_overlaps.primer_set_id,
       variant_overlaps.primer_set_name,
       variant_overlaps.variant_id,
       variant_overlaps.variant_type,
       variant_overlaps.variant,
       variant_overlaps.variant_start,
       variant_overlaps.variant_end,
       variant_overlaps.region,
       variant_overlaps.subregion,
       variant_overlaps.division,
       variant_overlaps.subdivision,
       variant_overlaps.detailed_geo_location_id,
       variant_overlaps.date_collected,
       variant_overlaps.variant_start AS coords
FROM public.variant_overlaps
WHERE ((variant_overlaps.variant_type)::text = 'I'::text)
WITH NO DATA
;

CREATE INDEX ON oligo_variant_overlaps (oligo_id);
CREATE INDEX ON oligo_variant_overlaps (date_collected);
CREATE INDEX ON oligo_variant_overlaps (detailed_geo_location_id);
GRANT SELECT on oligo_variant_overlaps to primer_monitor_ro;

ALTER MATERIALIZED VIEW oligo_variant_overlaps owner to primer_monitor;



