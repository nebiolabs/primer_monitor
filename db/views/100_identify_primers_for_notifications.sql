

CREATE MATERIALIZED VIEW public.identify_primers_for_notifications AS
 WITH first_query AS (
         SELECT primer_set_subscriptions.user_id,
            primer_set_subscriptions.primer_set_id,
            join_subscribed_location_to_ids.detailed_geo_location_id,
            primer_sets.name AS set_name,
            oligos.name AS primer_name,
            oligos.id AS oligo_id,
            users.lookback_days,
            users.variant_fraction_threshold,
            oligo_variant_overlaps.region,
            oligo_variant_overlaps.subregion,
            oligo_variant_overlaps.division,
            oligo_variant_overlaps.subdivision,
            oligo_variant_overlaps.detailed_geo_location_id AS unused_id,
            oligo_variant_overlaps.coords,
            count(oligo_variant_overlaps.variant_id) AS variant_count
           FROM (((((public.primer_set_subscriptions
             JOIN public.primer_sets ON ((primer_sets.id = primer_set_subscriptions.primer_set_id)))
             JOIN public.oligos ON ((primer_sets.id = oligos.primer_set_id)))
             JOIN public.oligo_variant_overlaps ON ((oligo_variant_overlaps.oligo_id = oligos.id)))
             JOIN public.join_subscribed_location_to_ids ON (((join_subscribed_location_to_ids.user_id = primer_set_subscriptions.user_id) AND (join_subscribed_location_to_ids.detailed_geo_location_id = oligo_variant_overlaps.detailed_geo_location_id))))
             JOIN public.users ON ((users.id = primer_set_subscriptions.user_id)))
          WHERE (oligo_variant_overlaps.date_collected >= (CURRENT_DATE - users.lookback_days))
          GROUP BY primer_set_subscriptions.user_id, primer_set_subscriptions.primer_set_id, primer_sets.name, oligos.id, oligos.name, join_subscribed_location_to_ids.detailed_geo_location_id, users.lookback_days, users.variant_fraction_threshold, oligo_variant_overlaps.region, oligo_variant_overlaps.subregion, oligo_variant_overlaps.division, oligo_variant_overlaps.subdivision, oligo_variant_overlaps.coords, oligo_variant_overlaps.detailed_geo_location_id
        ), second_query AS (
         SELECT fasta_records.detailed_geo_location_id,
            count(fasta_records.id) AS records_count,
            users.lookback_days
           FROM ((public.fasta_records
             JOIN public.join_subscribed_location_to_ids ON ((join_subscribed_location_to_ids.detailed_geo_location_id = fasta_records.detailed_geo_location_id)))
             JOIN public.users ON ((join_subscribed_location_to_ids.user_id = users.id)))
          WHERE (fasta_records.date_collected >= (CURRENT_DATE - users.lookback_days))
          GROUP BY fasta_records.detailed_geo_location_id, users.lookback_days
         HAVING (count(fasta_records.id) >= 20)
        )
 SELECT first_query.user_id,
    first_query.primer_set_id,
    first_query.set_name,
    first_query.oligo_id,
    first_query.primer_name,
    first_query.region,
    first_query.subregion,
    first_query.division,
    first_query.subdivision,
    first_query.coords,
    first_query.variant_count,
    first_query.variant_fraction_threshold,
    second_query.detailed_geo_location_id,
    second_query.records_count,
    ((first_query.variant_count)::numeric / (second_query.records_count)::numeric) AS fraction_variant
   FROM (first_query
     JOIN second_query ON ((second_query.detailed_geo_location_id = first_query.detailed_geo_location_id)))
  WHERE ((((first_query.variant_count)::numeric / (second_query.records_count)::numeric))::double precision >= first_query.variant_fraction_threshold)
  ;

