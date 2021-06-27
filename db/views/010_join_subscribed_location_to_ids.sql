
CREATE VIEW public.join_subscribed_location_to_ids AS
 WITH subscribed_ids AS (
         SELECT subscribed_geo_locations.user_id,
            subscribed_geo_locations.detailed_geo_location_alias_id,
            detailed_geo_location_aliases.region,
            detailed_geo_location_aliases.subregion,
            detailed_geo_location_aliases.division,
            detailed_geo_location_aliases.subdivision
           FROM (public.detailed_geo_location_aliases
             JOIN public.subscribed_geo_locations ON ((subscribed_geo_locations.detailed_geo_location_alias_id = detailed_geo_location_aliases.id)))
        )
 SELECT subscribed_ids.user_id,
    detailed_geo_locations.id AS detailed_geo_location_id,
    subscribed_ids.detailed_geo_location_alias_id
   FROM (public.detailed_geo_locations
     JOIN subscribed_ids ON ((((subscribed_ids.region IS NULL) OR ((subscribed_ids.region)::text = (detailed_geo_locations.region)::text)) AND ((subscribed_ids.subregion IS NULL) OR ((subscribed_ids.subregion)::text = (detailed_geo_locations.subregion)::text)) AND ((subscribed_ids.division IS NULL) OR ((subscribed_ids.division)::text = (detailed_geo_locations.division)::text)) AND ((subscribed_ids.subdivision IS NULL) OR ((subscribed_ids.subdivision)::text = (detailed_geo_locations.subdivision)::text)))));
;

