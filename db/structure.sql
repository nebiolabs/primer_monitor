SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET xmloption = content;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: oligo_category; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.oligo_category AS ENUM (
    'BIP',
    'FIP',
    'LB',
    'LF',
    'B3',
    'F3',
    'Reverse',
    'Forward',
    'Probe'
);


--
-- Name: primer_set_status; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE public.primer_set_status AS ENUM (
    'pending',
    'complete',
    'failed'
);


SET default_tablespace = '';

SET default_table_access_method = heap;

--
-- Name: amplification_methods; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.amplification_methods (
    id bigint NOT NULL,
    name character varying,
    description_url character varying,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
);


--
-- Name: amplification_methods_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.amplification_methods_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: amplification_methods_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.amplification_methods_id_seq OWNED BY public.amplification_methods.id;


--
-- Name: ar_internal_metadata; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.ar_internal_metadata (
    key character varying NOT NULL,
    value character varying,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
);


--
-- Name: blast_hits; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.blast_hits (
    id bigint NOT NULL,
    oligo_id bigint,
    organism_id bigint,
    num_identities integer,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
);


--
-- Name: blast_hits_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.blast_hits_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: blast_hits_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.blast_hits_id_seq OWNED BY public.blast_hits.id;


--
-- Name: detailed_geo_location_aliases; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.detailed_geo_location_aliases (
    id bigint NOT NULL,
    world character varying,
    region character varying,
    subregion character varying,
    division character varying,
    subdivision character varying,
    locality character varying,
    sublocality character varying,
    latitude double precision,
    longitude double precision,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
);


--
-- Name: detailed_geo_location_aliases_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.detailed_geo_location_aliases_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: detailed_geo_location_aliases_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.detailed_geo_location_aliases_id_seq OWNED BY public.detailed_geo_location_aliases.id;


--
-- Name: detailed_geo_locations; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.detailed_geo_locations (
    id bigint NOT NULL,
    world character varying NOT NULL,
    region character varying,
    subregion character varying,
    division character varying,
    subdivision character varying,
    locality character varying,
    sublocality character varying,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
);


--
-- Name: detailed_geo_locations_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.detailed_geo_locations_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: detailed_geo_locations_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.detailed_geo_locations_id_seq OWNED BY public.detailed_geo_locations.id;


--
-- Name: fasta_records; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.fasta_records (
    id bigint NOT NULL,
    strain character varying,
    genbank_accession character varying,
    gisaid_epi_isl character varying,
    date_collected date,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    variant_name character varying,
    detailed_geo_location_id bigint
);


--
-- Name: fasta_records_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.fasta_records_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: fasta_records_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.fasta_records_id_seq OWNED BY public.fasta_records.id;


--
-- Name: genomic_features; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.genomic_features (
    id bigint NOT NULL,
    name character varying,
    type character varying,
    ref_start integer,
    ref_end integer,
    organism_id bigint NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
);


--
-- Name: genomic_features_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.genomic_features_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: genomic_features_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.genomic_features_id_seq OWNED BY public.genomic_features.id;


--
-- Name: subscribed_geo_locations; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.subscribed_geo_locations (
    id bigint NOT NULL,
    user_id bigint NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    detailed_geo_location_alias_id bigint
);


--
-- Name: join_subscribed_location_to_id; Type: VIEW; Schema: public; Owner: -
--

CREATE VIEW public.join_subscribed_location_to_id AS
 WITH subscribed_ids AS (
         SELECT subscribed_geo_locations.user_id,
            detailed_geo_location_aliases.region,
            detailed_geo_location_aliases.subregion,
            detailed_geo_location_aliases.division,
            detailed_geo_location_aliases.subdivision
           FROM (public.detailed_geo_location_aliases
             JOIN public.subscribed_geo_locations ON ((subscribed_geo_locations.detailed_geo_location_alias_id = detailed_geo_location_aliases.id)))
        )
 SELECT subscribed_ids.user_id,
    detailed_geo_locations.id AS detailed_geo_location_id
   FROM (public.detailed_geo_locations
     JOIN subscribed_ids ON ((((subscribed_ids.region IS NULL) OR ((subscribed_ids.region)::text = (detailed_geo_locations.region)::text)) AND ((subscribed_ids.subregion IS NULL) OR ((subscribed_ids.subregion)::text = (detailed_geo_locations.subregion)::text)) AND ((subscribed_ids.division IS NULL) OR ((subscribed_ids.division)::text = (detailed_geo_locations.division)::text)) AND ((subscribed_ids.subdivision IS NULL) OR ((subscribed_ids.subdivision)::text = (detailed_geo_locations.subdivision)::text)))));


--
-- Name: oligos; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.oligos (
    id bigint NOT NULL,
    name character varying NOT NULL,
    sequence character varying NOT NULL,
    primer_set_id bigint NOT NULL,
    ref_start bigint,
    ref_end bigint,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    locus character varying,
    category public.oligo_category,
    short_name character varying
);


--
-- Name: variant_sites; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.variant_sites (
    id bigint NOT NULL,
    ref_start integer,
    variant_type character varying,
    variant character varying,
    fasta_record_id bigint,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    ref_end integer
);


--
-- Name: oligo_variant_overlaps; Type: MATERIALIZED VIEW; Schema: public; Owner: -
--

CREATE MATERIALIZED VIEW public.oligo_variant_overlaps AS
 WITH big_query AS (
         SELECT oligos.id AS oligo_id,
            oligos.name AS oligo_name,
            oligos.ref_start AS oligo_start,
            oligos.ref_end AS oligo_end,
            variant_sites.id AS variant_id,
            variant_sites.variant_type,
            variant_sites.variant,
            variant_sites.ref_start AS variant_start,
            variant_sites.ref_end AS variant_end,
            detailed_geo_locations.region,
            detailed_geo_locations.subregion,
            detailed_geo_locations.division,
            detailed_geo_locations.subdivision,
            detailed_geo_locations.id AS detailed_geo_location_id,
            fasta_records.date_collected
           FROM (((public.variant_sites
             JOIN public.fasta_records ON ((variant_sites.fasta_record_id = fasta_records.id)))
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
             JOIN public.oligos ON ((((variant_sites.ref_start >= oligos.ref_start) AND (variant_sites.ref_start <= oligos.ref_end)) OR ((variant_sites.ref_end >= oligos.ref_start) AND (variant_sites.ref_end <= oligos.ref_end)) OR ((variant_sites.ref_start < oligos.ref_start) AND (variant_sites.ref_end > oligos.ref_end)))))
          WHERE ((((variant_sites.variant_type)::text = 'D'::text) OR ((variant_sites.variant_type)::text = 'X'::text)) AND ((variant_sites.variant)::text !~~ '%N%'::text))
        ), insert_query AS (
         SELECT oligos.id AS oligo_id,
            oligos.name AS oligo_name,
            oligos.ref_start AS oligo_start,
            oligos.ref_end AS oligo_end,
            variant_sites.id AS variant_id,
            variant_sites.variant_type,
            variant_sites.variant,
            variant_sites.ref_start AS variant_start,
            variant_sites.ref_end AS variant_end,
            detailed_geo_locations.region,
            detailed_geo_locations.subregion,
            detailed_geo_locations.division,
            detailed_geo_locations.subdivision,
            detailed_geo_locations.id AS detailed_geo_location_id,
            fasta_records.date_collected
           FROM (((public.variant_sites
             JOIN public.fasta_records ON ((variant_sites.fasta_record_id = fasta_records.id)))
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
             JOIN public.oligos ON ((((variant_sites.ref_start >= oligos.ref_start) AND (variant_sites.ref_start <= oligos.ref_end)) OR ((variant_sites.ref_end >= oligos.ref_start) AND (variant_sites.ref_end <= oligos.ref_end)) OR ((variant_sites.ref_start < oligos.ref_start) AND (variant_sites.ref_end > oligos.ref_end)))))
          WHERE (((variant_sites.variant_type)::text = 'I'::text) AND ((variant_sites.variant)::text !~~ '%N%'::text))
        ), region_count AS (
         SELECT count(*) AS region_count,
            detailed_geo_locations.region
           FROM (public.fasta_records
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
          GROUP BY detailed_geo_locations.region
        ), region_subregion_count AS (
         SELECT count(*) AS region_subregion_count,
            detailed_geo_locations.region,
            detailed_geo_locations.subregion
           FROM (public.fasta_records
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion
        ), region_subregion_division_count AS (
         SELECT count(*) AS region_subregion_division_count,
            detailed_geo_locations.region,
            detailed_geo_locations.subregion,
            detailed_geo_locations.division
           FROM (public.fasta_records
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division
        ), region_subregion_division_subdivision_count AS (
         SELECT count(*) AS region_subregion_division_subdivision_count,
            detailed_geo_locations.region,
            detailed_geo_locations.subregion,
            detailed_geo_locations.division,
            detailed_geo_locations.subdivision
           FROM (public.fasta_records
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, detailed_geo_locations.subdivision
        ), region_time_count AS (
         SELECT count(*) AS region_time_count,
            detailed_geo_locations.region,
            fasta_records.date_collected
           FROM (public.fasta_records
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
          GROUP BY detailed_geo_locations.region, fasta_records.date_collected
        ), region_subregion_time_count AS (
         SELECT count(*) AS region_subregion_time_count,
            detailed_geo_locations.region,
            detailed_geo_locations.subregion,
            fasta_records.date_collected
           FROM (public.fasta_records
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, fasta_records.date_collected
        ), region_subregion_division_time_count AS (
         SELECT count(*) AS region_subregion_division_time_count,
            detailed_geo_locations.region,
            detailed_geo_locations.subregion,
            detailed_geo_locations.division,
            fasta_records.date_collected
           FROM (public.fasta_records
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, fasta_records.date_collected
        ), region_subregion_division_subdivision_time_count AS (
         SELECT count(*) AS region_subregion_division_subdivision_time_count,
            detailed_geo_locations.region,
            detailed_geo_locations.subregion,
            detailed_geo_locations.division,
            detailed_geo_locations.subdivision,
            fasta_records.date_collected
           FROM (public.fasta_records
             JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
          GROUP BY detailed_geo_locations.region, detailed_geo_locations.subregion, detailed_geo_locations.division, detailed_geo_locations.subdivision, fasta_records.date_collected
        )
 SELECT big_query.oligo_id,
    big_query.oligo_name,
    big_query.oligo_start,
    big_query.oligo_end,
    big_query.variant_id,
    big_query.variant_type,
    big_query.variant,
    big_query.variant_start,
    big_query.variant_end,
    big_query.region,
    big_query.subregion,
    big_query.division,
    big_query.subdivision,
    big_query.detailed_geo_location_id,
    big_query.date_collected,
    region_count.region_count,
    region_subregion_count.region_subregion_count,
    region_subregion_division_count.region_subregion_division_count,
    region_subregion_division_subdivision_count.region_subregion_division_subdivision_count,
    region_time_count.region_time_count,
    region_subregion_time_count.region_subregion_time_count,
    region_subregion_division_time_count.region_subregion_division_time_count,
    region_subregion_division_subdivision_time_count.region_subregion_division_subdivision_time_count,
    generate_series(lower((numrange((coord_overlaps.oligo_start)::numeric, (coord_overlaps.oligo_end)::numeric) * numrange((coord_overlaps.variant_start)::numeric, (coord_overlaps.variant_end)::numeric))), (upper((numrange((coord_overlaps.oligo_start)::numeric, (coord_overlaps.oligo_end)::numeric) * numrange((coord_overlaps.variant_start)::numeric, (coord_overlaps.variant_end)::numeric))) - (1)::numeric)) AS coords
   FROM (((((((((big_query coord_overlaps
     JOIN big_query ON (((coord_overlaps.oligo_id = big_query.oligo_id) AND (coord_overlaps.variant_id = big_query.variant_id))))
     JOIN region_count ON (((region_count.region)::text = (big_query.region)::text)))
     JOIN region_subregion_count ON ((((region_subregion_count.region)::text = (big_query.region)::text) AND ((region_subregion_count.subregion)::text = (big_query.subregion)::text))))
     JOIN region_subregion_division_count ON ((((region_subregion_division_count.region)::text = (big_query.region)::text) AND ((region_subregion_division_count.subregion)::text = (big_query.subregion)::text) AND ((region_subregion_division_count.division)::text = (big_query.division)::text))))
     JOIN region_subregion_division_subdivision_count ON ((((region_subregion_division_subdivision_count.region)::text = (big_query.region)::text) AND ((region_subregion_division_subdivision_count.subregion)::text = (big_query.subregion)::text) AND ((region_subregion_division_subdivision_count.division)::text = (big_query.division)::text) AND ((region_subregion_division_subdivision_count.subdivision)::text = (big_query.subdivision)::text))))
     JOIN region_time_count ON ((((region_time_count.region)::text = (big_query.region)::text) AND (region_time_count.date_collected = big_query.date_collected))))
     JOIN region_subregion_time_count ON ((((region_subregion_time_count.region)::text = (big_query.region)::text) AND ((region_subregion_time_count.subregion)::text = (big_query.subregion)::text) AND (region_subregion_time_count.date_collected = big_query.date_collected))))
     JOIN region_subregion_division_time_count ON ((((region_subregion_division_time_count.region)::text = (big_query.region)::text) AND ((region_subregion_division_time_count.subregion)::text = (big_query.subregion)::text) AND ((region_subregion_division_time_count.division)::text = (big_query.division)::text) AND (region_subregion_division_time_count.date_collected = big_query.date_collected))))
     JOIN region_subregion_division_subdivision_time_count ON ((((region_subregion_division_subdivision_time_count.region)::text = (big_query.region)::text) AND ((region_subregion_division_subdivision_time_count.subregion)::text = (big_query.subregion)::text) AND ((region_subregion_division_subdivision_time_count.division)::text = (big_query.division)::text) AND ((region_subregion_division_subdivision_time_count.subdivision)::text = (big_query.subdivision)::text) AND (region_subregion_division_subdivision_time_count.date_collected = big_query.date_collected))))
UNION ALL
 SELECT insert_query.oligo_id,
    insert_query.oligo_name,
    insert_query.oligo_start,
    insert_query.oligo_end,
    insert_query.variant_id,
    insert_query.variant_type,
    insert_query.variant,
    insert_query.variant_start,
    insert_query.variant_end,
    insert_query.region,
    insert_query.subregion,
    insert_query.division,
    insert_query.subdivision,
    insert_query.detailed_geo_location_id,
    insert_query.date_collected,
    region_count.region_count,
    region_subregion_count.region_subregion_count,
    region_subregion_division_count.region_subregion_division_count,
    region_subregion_division_subdivision_count.region_subregion_division_subdivision_count,
    region_time_count.region_time_count,
    region_subregion_time_count.region_subregion_time_count,
    region_subregion_division_time_count.region_subregion_division_time_count,
    region_subregion_division_subdivision_time_count.region_subregion_division_subdivision_time_count,
    insert_query.variant_start AS coords
   FROM (((((((((insert_query coord_overlaps
     JOIN insert_query ON (((coord_overlaps.oligo_id = insert_query.oligo_id) AND (coord_overlaps.variant_id = insert_query.variant_id))))
     JOIN region_count ON (((region_count.region)::text = (insert_query.region)::text)))
     JOIN region_subregion_count ON ((((region_subregion_count.region)::text = (insert_query.region)::text) AND ((region_subregion_count.subregion)::text = (insert_query.subregion)::text))))
     JOIN region_subregion_division_count ON ((((region_subregion_division_count.region)::text = (insert_query.region)::text) AND ((region_subregion_division_count.subregion)::text = (insert_query.subregion)::text) AND ((region_subregion_division_count.division)::text = (insert_query.division)::text))))
     JOIN region_subregion_division_subdivision_count ON ((((region_subregion_division_subdivision_count.region)::text = (insert_query.region)::text) AND ((region_subregion_division_subdivision_count.subregion)::text = (insert_query.subregion)::text) AND ((region_subregion_division_subdivision_count.division)::text = (insert_query.division)::text) AND ((region_subregion_division_subdivision_count.subdivision)::text = (insert_query.subdivision)::text))))
     JOIN region_time_count ON ((((region_time_count.region)::text = (insert_query.region)::text) AND (region_time_count.date_collected = insert_query.date_collected))))
     JOIN region_subregion_time_count ON ((((region_subregion_time_count.region)::text = (insert_query.region)::text) AND ((region_subregion_time_count.subregion)::text = (insert_query.subregion)::text) AND (region_subregion_time_count.date_collected = insert_query.date_collected))))
     JOIN region_subregion_division_time_count ON ((((region_subregion_division_time_count.region)::text = (insert_query.region)::text) AND ((region_subregion_division_time_count.subregion)::text = (insert_query.subregion)::text) AND ((region_subregion_division_time_count.division)::text = (insert_query.division)::text) AND (region_subregion_division_time_count.date_collected = insert_query.date_collected))))
     JOIN region_subregion_division_subdivision_time_count ON ((((region_subregion_division_subdivision_time_count.region)::text = (insert_query.region)::text) AND ((region_subregion_division_subdivision_time_count.subregion)::text = (insert_query.subregion)::text) AND ((region_subregion_division_subdivision_time_count.division)::text = (insert_query.division)::text) AND ((region_subregion_division_subdivision_time_count.subdivision)::text = (insert_query.subdivision)::text) AND (region_subregion_division_subdivision_time_count.date_collected = insert_query.date_collected))))
  WITH NO DATA;


--
-- Name: primer_set_subscriptions; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.primer_set_subscriptions (
    id bigint NOT NULL,
    primer_set_id bigint,
    user_id bigint,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
);


--
-- Name: primer_sets; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.primer_sets (
    id bigint NOT NULL,
    name character varying,
    user_id bigint,
    organism_id integer NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    status public.primer_set_status DEFAULT 'pending'::public.primer_set_status,
    citation_url character varying,
    doi character varying,
    amplification_method_id bigint
);


--
-- Name: users; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.users (
    id bigint NOT NULL,
    first character varying NOT NULL,
    last character varying NOT NULL,
    email character varying NOT NULL,
    activated boolean DEFAULT false NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    login character varying,
    encrypted_password character varying,
    sign_in_count integer DEFAULT 0 NOT NULL,
    failed_login_count integer DEFAULT 0 NOT NULL,
    last_request_at timestamp without time zone,
    current_sign_in_at timestamp without time zone,
    last_sign_in_at timestamp without time zone,
    current_sign_in_ip character varying,
    last_sign_in_ip character varying,
    active boolean DEFAULT false,
    approved boolean DEFAULT false,
    confirmed boolean DEFAULT false,
    send_primer_updates boolean DEFAULT false NOT NULL,
    reset_password_token character varying,
    reset_password_sent_at timestamp without time zone,
    remember_token character varying,
    remember_created_at timestamp without time zone,
    authentication_token character varying,
    confirmation_token character varying(255),
    confirmed_at timestamp without time zone,
    confirmation_sent_at timestamp without time zone,
    unconfirmed_email character varying,
    failed_attempts integer DEFAULT 0 NOT NULL,
    unlock_token character varying,
    locked_at timestamp without time zone,
    provider character varying,
    uid character varying,
    lookback_days integer DEFAULT 30 NOT NULL,
    variant_fraction_threshold double precision DEFAULT 0.1 NOT NULL
);


--
-- Name: identify_primers_for_notification; Type: MATERIALIZED VIEW; Schema: public; Owner: -
--

CREATE MATERIALIZED VIEW public.identify_primers_for_notification AS
 WITH first_query AS (
         SELECT primer_set_subscriptions.user_id,
            primer_set_subscriptions.primer_set_id,
            join_subscribed_location_to_id.detailed_geo_location_id,
            primer_sets.name AS set_name,
            oligos.name AS primer_name,
            users.lookback_days,
            users.variant_fraction_threshold,
            oligo_variant_overlaps.region,
            oligo_variant_overlaps.division,
            oligo_variant_overlaps.detailed_geo_location_id AS unused_id,
            oligo_variant_overlaps.coords,
            count(oligo_variant_overlaps.variant_id) AS variant_count
           FROM (((((public.primer_set_subscriptions
             JOIN public.primer_sets ON ((primer_sets.id = primer_set_subscriptions.primer_set_id)))
             JOIN public.join_subscribed_location_to_id ON ((join_subscribed_location_to_id.user_id = primer_set_subscriptions.user_id)))
             JOIN public.oligos ON ((primer_sets.id = oligos.primer_set_id)))
             JOIN public.oligo_variant_overlaps ON ((oligo_variant_overlaps.oligo_id = oligos.id)))
             JOIN public.users ON ((users.id = primer_set_subscriptions.user_id)))
          WHERE (oligo_variant_overlaps.date_collected >= (CURRENT_DATE - users.lookback_days))
          GROUP BY primer_set_subscriptions.user_id, primer_set_subscriptions.primer_set_id, primer_sets.name, oligos.name, join_subscribed_location_to_id.detailed_geo_location_id, users.lookback_days, users.variant_fraction_threshold, oligo_variant_overlaps.region, oligo_variant_overlaps.division, oligo_variant_overlaps.coords, oligo_variant_overlaps.detailed_geo_location_id
        ), second_query AS (
         SELECT fasta_records.detailed_geo_location_id,
            count(fasta_records.id) AS records_count,
            users.lookback_days
           FROM ((public.fasta_records
             JOIN public.join_subscribed_location_to_id ON ((join_subscribed_location_to_id.detailed_geo_location_id = fasta_records.detailed_geo_location_id)))
             JOIN public.users ON ((join_subscribed_location_to_id.user_id = users.id)))
          WHERE (fasta_records.date_collected >= (CURRENT_DATE - users.lookback_days))
          GROUP BY fasta_records.detailed_geo_location_id, users.lookback_days
         HAVING (count(fasta_records.id) >= 20)
        )
 SELECT first_query.user_id,
    first_query.primer_set_id,
    first_query.set_name,
    first_query.primer_name,
    first_query.region,
    first_query.division,
    first_query.coords,
    first_query.variant_count,
    first_query.variant_fraction_threshold,
    second_query.detailed_geo_location_id,
    second_query.records_count,
    ((first_query.variant_count)::numeric / (second_query.records_count)::numeric) AS fraction_variant
   FROM (first_query
     JOIN second_query ON ((second_query.detailed_geo_location_id = first_query.detailed_geo_location_id)))
  WHERE ((((first_query.variant_count)::numeric / (second_query.records_count)::numeric))::double precision >= first_query.variant_fraction_threshold)
  WITH NO DATA;


--
-- Name: location_alias_join; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.location_alias_join (
    id bigint NOT NULL,
    detailed_geo_locations_id bigint,
    detailed_geo_location_aliases_id bigint,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
);


--
-- Name: location_alias_join_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.location_alias_join_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: location_alias_join_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.location_alias_join_id_seq OWNED BY public.location_alias_join.id;


--
-- Name: oligos_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.oligos_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: oligos_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.oligos_id_seq OWNED BY public.oligos.id;


--
-- Name: organisms; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.organisms (
    id bigint NOT NULL,
    ncbi_taxon_id integer NOT NULL,
    name character varying NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    alias character varying
);


--
-- Name: organisms_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.organisms_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: organisms_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.organisms_id_seq OWNED BY public.organisms.id;


--
-- Name: primer_set_subscriptions_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.primer_set_subscriptions_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: primer_set_subscriptions_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.primer_set_subscriptions_id_seq OWNED BY public.primer_set_subscriptions.id;


--
-- Name: primer_sets_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.primer_sets_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: primer_sets_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.primer_sets_id_seq OWNED BY public.primer_sets.id;


--
-- Name: proposed_notifications; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.proposed_notifications (
    id bigint NOT NULL,
    primer_sets_id bigint,
    oligos_id bigint,
    verified_notifications_id bigint,
    coordinate integer,
    fraction_variant double precision,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
);


--
-- Name: proposed_notifications_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.proposed_notifications_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: proposed_notifications_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.proposed_notifications_id_seq OWNED BY public.proposed_notifications.id;


--
-- Name: roles; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.roles (
    id bigint NOT NULL,
    name character varying,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
);


--
-- Name: roles_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.roles_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: roles_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.roles_id_seq OWNED BY public.roles.id;


--
-- Name: schema_migrations; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.schema_migrations (
    version character varying NOT NULL
);


--
-- Name: subscribed_geo_locations_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.subscribed_geo_locations_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: subscribed_geo_locations_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.subscribed_geo_locations_id_seq OWNED BY public.subscribed_geo_locations.id;


--
-- Name: user_roles; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.user_roles (
    id bigint NOT NULL,
    user_id bigint NOT NULL,
    role_id bigint NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
);


--
-- Name: user_roles_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.user_roles_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: user_roles_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.user_roles_id_seq OWNED BY public.user_roles.id;


--
-- Name: users_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.users_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: users_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.users_id_seq OWNED BY public.users.id;


--
-- Name: variant_sites_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.variant_sites_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: variant_sites_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.variant_sites_id_seq OWNED BY public.variant_sites.id;


--
-- Name: verified_notifications; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.verified_notifications (
    id bigint NOT NULL,
    users_id bigint,
    date_sent date,
    status character varying,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
);


--
-- Name: verified_notifications_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.verified_notifications_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: verified_notifications_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.verified_notifications_id_seq OWNED BY public.verified_notifications.id;


--
-- Name: amplification_methods id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.amplification_methods ALTER COLUMN id SET DEFAULT nextval('public.amplification_methods_id_seq'::regclass);


--
-- Name: blast_hits id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.blast_hits ALTER COLUMN id SET DEFAULT nextval('public.blast_hits_id_seq'::regclass);


--
-- Name: detailed_geo_location_aliases id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.detailed_geo_location_aliases ALTER COLUMN id SET DEFAULT nextval('public.detailed_geo_location_aliases_id_seq'::regclass);


--
-- Name: detailed_geo_locations id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.detailed_geo_locations ALTER COLUMN id SET DEFAULT nextval('public.detailed_geo_locations_id_seq'::regclass);


--
-- Name: fasta_records id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.fasta_records ALTER COLUMN id SET DEFAULT nextval('public.fasta_records_id_seq'::regclass);


--
-- Name: genomic_features id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genomic_features ALTER COLUMN id SET DEFAULT nextval('public.genomic_features_id_seq'::regclass);


--
-- Name: location_alias_join id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.location_alias_join ALTER COLUMN id SET DEFAULT nextval('public.location_alias_join_id_seq'::regclass);


--
-- Name: oligos id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.oligos ALTER COLUMN id SET DEFAULT nextval('public.oligos_id_seq'::regclass);


--
-- Name: organisms id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.organisms ALTER COLUMN id SET DEFAULT nextval('public.organisms_id_seq'::regclass);


--
-- Name: primer_set_subscriptions id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.primer_set_subscriptions ALTER COLUMN id SET DEFAULT nextval('public.primer_set_subscriptions_id_seq'::regclass);


--
-- Name: primer_sets id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.primer_sets ALTER COLUMN id SET DEFAULT nextval('public.primer_sets_id_seq'::regclass);


--
-- Name: proposed_notifications id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.proposed_notifications ALTER COLUMN id SET DEFAULT nextval('public.proposed_notifications_id_seq'::regclass);


--
-- Name: roles id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.roles ALTER COLUMN id SET DEFAULT nextval('public.roles_id_seq'::regclass);


--
-- Name: subscribed_geo_locations id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.subscribed_geo_locations ALTER COLUMN id SET DEFAULT nextval('public.subscribed_geo_locations_id_seq'::regclass);


--
-- Name: user_roles id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.user_roles ALTER COLUMN id SET DEFAULT nextval('public.user_roles_id_seq'::regclass);


--
-- Name: users id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.users ALTER COLUMN id SET DEFAULT nextval('public.users_id_seq'::regclass);


--
-- Name: variant_sites id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variant_sites ALTER COLUMN id SET DEFAULT nextval('public.variant_sites_id_seq'::regclass);


--
-- Name: verified_notifications id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.verified_notifications ALTER COLUMN id SET DEFAULT nextval('public.verified_notifications_id_seq'::regclass);


--
-- Name: amplification_methods amplification_methods_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.amplification_methods
    ADD CONSTRAINT amplification_methods_pkey PRIMARY KEY (id);


--
-- Name: ar_internal_metadata ar_internal_metadata_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.ar_internal_metadata
    ADD CONSTRAINT ar_internal_metadata_pkey PRIMARY KEY (key);


--
-- Name: blast_hits blast_hits_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.blast_hits
    ADD CONSTRAINT blast_hits_pkey PRIMARY KEY (id);


--
-- Name: detailed_geo_location_aliases detailed_geo_location_aliases_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.detailed_geo_location_aliases
    ADD CONSTRAINT detailed_geo_location_aliases_pkey PRIMARY KEY (id);


--
-- Name: detailed_geo_locations detailed_geo_locations_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.detailed_geo_locations
    ADD CONSTRAINT detailed_geo_locations_pkey PRIMARY KEY (id);


--
-- Name: fasta_records fasta_records_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.fasta_records
    ADD CONSTRAINT fasta_records_pkey PRIMARY KEY (id);


--
-- Name: genomic_features genomic_features_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genomic_features
    ADD CONSTRAINT genomic_features_pkey PRIMARY KEY (id);


--
-- Name: location_alias_join location_alias_join_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.location_alias_join
    ADD CONSTRAINT location_alias_join_pkey PRIMARY KEY (id);


--
-- Name: oligos oligos_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.oligos
    ADD CONSTRAINT oligos_pkey PRIMARY KEY (id);


--
-- Name: organisms organisms_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.organisms
    ADD CONSTRAINT organisms_pkey PRIMARY KEY (id);


--
-- Name: primer_set_subscriptions primer_set_subscriptions_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.primer_set_subscriptions
    ADD CONSTRAINT primer_set_subscriptions_pkey PRIMARY KEY (id);


--
-- Name: primer_sets primer_sets_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.primer_sets
    ADD CONSTRAINT primer_sets_pkey PRIMARY KEY (id);


--
-- Name: proposed_notifications proposed_notifications_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.proposed_notifications
    ADD CONSTRAINT proposed_notifications_pkey PRIMARY KEY (id);


--
-- Name: roles roles_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.roles
    ADD CONSTRAINT roles_pkey PRIMARY KEY (id);


--
-- Name: schema_migrations schema_migrations_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.schema_migrations
    ADD CONSTRAINT schema_migrations_pkey PRIMARY KEY (version);


--
-- Name: subscribed_geo_locations subscribed_geo_locations_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.subscribed_geo_locations
    ADD CONSTRAINT subscribed_geo_locations_pkey PRIMARY KEY (id);


--
-- Name: user_roles user_roles_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.user_roles
    ADD CONSTRAINT user_roles_pkey PRIMARY KEY (id);


--
-- Name: users users_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.users
    ADD CONSTRAINT users_pkey PRIMARY KEY (id);


--
-- Name: variant_sites variant_sites_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variant_sites
    ADD CONSTRAINT variant_sites_pkey PRIMARY KEY (id);


--
-- Name: verified_notifications verified_notifications_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.verified_notifications
    ADD CONSTRAINT verified_notifications_pkey PRIMARY KEY (id);


--
-- Name: alias_full_record; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX alias_full_record ON public.detailed_geo_location_aliases USING btree (region, subregion, division, subdivision, locality, sublocality, latitude, longitude);


--
-- Name: alias_join_index; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX alias_join_index ON public.location_alias_join USING btree (detailed_geo_locations_id);


--
-- Name: full_record; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX full_record ON public.detailed_geo_locations USING btree (region, subregion, division, subdivision, locality, sublocality);


--
-- Name: index_blast_hits_on_oligo_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_blast_hits_on_oligo_id ON public.blast_hits USING btree (oligo_id);


--
-- Name: index_blast_hits_on_organism_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_blast_hits_on_organism_id ON public.blast_hits USING btree (organism_id);


--
-- Name: index_fasta_records_on_detailed_geo_location_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_fasta_records_on_detailed_geo_location_id ON public.fasta_records USING btree (detailed_geo_location_id);


--
-- Name: index_fasta_records_on_strain; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX index_fasta_records_on_strain ON public.fasta_records USING btree (strain);


--
-- Name: index_genomic_features_on_organism_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_genomic_features_on_organism_id ON public.genomic_features USING btree (organism_id);


--
-- Name: index_location_alias_join_on_detailed_geo_locations_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_location_alias_join_on_detailed_geo_locations_id ON public.location_alias_join USING btree (detailed_geo_locations_id);


--
-- Name: index_oligos_on_primer_set_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_oligos_on_primer_set_id ON public.oligos USING btree (primer_set_id);


--
-- Name: index_primer_set_subscriptions_on_primer_set_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_primer_set_subscriptions_on_primer_set_id ON public.primer_set_subscriptions USING btree (primer_set_id);


--
-- Name: index_primer_set_subscriptions_on_user_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_primer_set_subscriptions_on_user_id ON public.primer_set_subscriptions USING btree (user_id);


--
-- Name: index_primer_set_subscriptions_on_user_id_and_primer_set_id; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX index_primer_set_subscriptions_on_user_id_and_primer_set_id ON public.primer_set_subscriptions USING btree (user_id, primer_set_id);


--
-- Name: index_primer_sets_on_organism_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_primer_sets_on_organism_id ON public.primer_sets USING btree (organism_id);


--
-- Name: index_primer_sets_on_user_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_primer_sets_on_user_id ON public.primer_sets USING btree (user_id);


--
-- Name: index_proposed_notifications_on_oligos_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_proposed_notifications_on_oligos_id ON public.proposed_notifications USING btree (oligos_id);


--
-- Name: index_proposed_notifications_on_primer_sets_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_proposed_notifications_on_primer_sets_id ON public.proposed_notifications USING btree (primer_sets_id);


--
-- Name: index_proposed_notifications_on_verified_notifications_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_proposed_notifications_on_verified_notifications_id ON public.proposed_notifications USING btree (verified_notifications_id);


--
-- Name: index_subscribed_geo_locations_on_user_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_subscribed_geo_locations_on_user_id ON public.subscribed_geo_locations USING btree (user_id);


--
-- Name: index_user_roles_on_role_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_user_roles_on_role_id ON public.user_roles USING btree (role_id);


--
-- Name: index_user_roles_on_user_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_user_roles_on_user_id ON public.user_roles USING btree (user_id);


--
-- Name: index_users_on_confirmation_token; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX index_users_on_confirmation_token ON public.users USING btree (confirmation_token);


--
-- Name: index_users_on_email; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX index_users_on_email ON public.users USING btree (email);


--
-- Name: index_users_on_reset_password_token; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX index_users_on_reset_password_token ON public.users USING btree (reset_password_token);


--
-- Name: index_users_on_unconfirmed_email; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX index_users_on_unconfirmed_email ON public.users USING btree (unconfirmed_email);


--
-- Name: index_users_on_unlock_token; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX index_users_on_unlock_token ON public.users USING btree (unlock_token);


--
-- Name: index_variant_sites_on_fasta_record_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_variant_sites_on_fasta_record_id ON public.variant_sites USING btree (fasta_record_id);


--
-- Name: index_verified_notifications_on_users_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_verified_notifications_on_users_id ON public.verified_notifications USING btree (users_id);


--
-- Name: tmp; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX tmp ON public.subscribed_geo_locations USING btree (detailed_geo_location_alias_id);


--
-- Name: primer_sets fk_rails_06b1f0d34e; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.primer_sets
    ADD CONSTRAINT fk_rails_06b1f0d34e FOREIGN KEY (amplification_method_id) REFERENCES public.amplification_methods(id);


--
-- Name: blast_hits fk_rails_1f04a34db0; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.blast_hits
    ADD CONSTRAINT fk_rails_1f04a34db0 FOREIGN KEY (organism_id) REFERENCES public.organisms(id);


--
-- Name: variant_sites fk_rails_2c76e30b83; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.variant_sites
    ADD CONSTRAINT fk_rails_2c76e30b83 FOREIGN KEY (fasta_record_id) REFERENCES public.fasta_records(id);


--
-- Name: user_roles fk_rails_318345354e; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.user_roles
    ADD CONSTRAINT fk_rails_318345354e FOREIGN KEY (user_id) REFERENCES public.users(id);


--
-- Name: proposed_notifications fk_rails_32a2a94275; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.proposed_notifications
    ADD CONSTRAINT fk_rails_32a2a94275 FOREIGN KEY (primer_sets_id) REFERENCES public.primer_sets(id);


--
-- Name: user_roles fk_rails_3369e0d5fc; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.user_roles
    ADD CONSTRAINT fk_rails_3369e0d5fc FOREIGN KEY (role_id) REFERENCES public.roles(id);


--
-- Name: proposed_notifications fk_rails_39f3a62e67; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.proposed_notifications
    ADD CONSTRAINT fk_rails_39f3a62e67 FOREIGN KEY (oligos_id) REFERENCES public.oligos(id);


--
-- Name: primer_sets fk_rails_52b9cf4012; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.primer_sets
    ADD CONSTRAINT fk_rails_52b9cf4012 FOREIGN KEY (organism_id) REFERENCES public.organisms(id);


--
-- Name: oligos fk_rails_58e7536518; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.oligos
    ADD CONSTRAINT fk_rails_58e7536518 FOREIGN KEY (primer_set_id) REFERENCES public.primer_sets(id);


--
-- Name: location_alias_join fk_rails_61222fa2c3; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.location_alias_join
    ADD CONSTRAINT fk_rails_61222fa2c3 FOREIGN KEY (detailed_geo_locations_id) REFERENCES public.detailed_geo_locations(id);


--
-- Name: genomic_features fk_rails_65e85371d4; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genomic_features
    ADD CONSTRAINT fk_rails_65e85371d4 FOREIGN KEY (organism_id) REFERENCES public.organisms(id);


--
-- Name: verified_notifications fk_rails_660c55653c; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.verified_notifications
    ADD CONSTRAINT fk_rails_660c55653c FOREIGN KEY (users_id) REFERENCES public.users(id);


--
-- Name: proposed_notifications fk_rails_736e9d2781; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.proposed_notifications
    ADD CONSTRAINT fk_rails_736e9d2781 FOREIGN KEY (verified_notifications_id) REFERENCES public.verified_notifications(id);


--
-- Name: subscribed_geo_locations fk_rails_7745c5f33b; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.subscribed_geo_locations
    ADD CONSTRAINT fk_rails_7745c5f33b FOREIGN KEY (user_id) REFERENCES public.users(id);


--
-- Name: subscribed_geo_locations fk_rails_7c2744b62d; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.subscribed_geo_locations
    ADD CONSTRAINT fk_rails_7c2744b62d FOREIGN KEY (detailed_geo_location_alias_id) REFERENCES public.detailed_geo_location_aliases(id);


--
-- Name: fasta_records fk_rails_82e2588d91; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.fasta_records
    ADD CONSTRAINT fk_rails_82e2588d91 FOREIGN KEY (detailed_geo_location_id) REFERENCES public.detailed_geo_locations(id);


--
-- Name: primer_set_subscriptions fk_rails_99041872f6; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.primer_set_subscriptions
    ADD CONSTRAINT fk_rails_99041872f6 FOREIGN KEY (primer_set_id) REFERENCES public.primer_sets(id);


--
-- Name: primer_sets fk_rails_a78cff2c70; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.primer_sets
    ADD CONSTRAINT fk_rails_a78cff2c70 FOREIGN KEY (user_id) REFERENCES public.users(id);


--
-- Name: blast_hits fk_rails_b63d58c7a5; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.blast_hits
    ADD CONSTRAINT fk_rails_b63d58c7a5 FOREIGN KEY (oligo_id) REFERENCES public.oligos(id);


--
-- Name: location_alias_join fk_rails_dff9f9948c; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.location_alias_join
    ADD CONSTRAINT fk_rails_dff9f9948c FOREIGN KEY (detailed_geo_location_aliases_id) REFERENCES public.detailed_geo_location_aliases(id);


--
-- Name: primer_set_subscriptions fk_rails_e7701775d5; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.primer_set_subscriptions
    ADD CONSTRAINT fk_rails_e7701775d5 FOREIGN KEY (user_id) REFERENCES public.users(id);


--
-- PostgreSQL database dump complete
--

SET search_path TO "$user", public;

INSERT INTO "schema_migrations" (version) VALUES
('20200809190433'),
('20200809190500'),
('20200809190541'),
('20200809190725'),
('20200810101557'),
('20201023152436'),
('20201023153637'),
('20201023161602'),
('20210125033513'),
('20210125033734'),
('20210125034331'),
('20210125101936'),
('20210128221802'),
('20210129031328'),
('20210130131727'),
('20210211193202'),
('20210211212104'),
('20210211220241'),
('20210212182903'),
('20210212195828'),
('20210212203816'),
('20210212204343'),
('20210216221118'),
('20210217145903'),
('20210217152635'),
('20210218005924'),
('20210218022343'),
('20210218123414'),
('20210218161052'),
('20210218161800'),
('20210218172439'),
('20210218172905'),
('20210218180015'),
('20210219000134'),
('20210219153745'),
('20210219185233'),
('20210220203358'),
('20210220215944'),
('20210220233440'),
('20210220234316'),
('20210220235805'),
('20210222122744'),
('20210223192700'),
('20210223194128'),
('20210224163709'),
('20210224192126'),
('20210225212556'),
('20210225225243'),
('20210226173653'),
('20210226174141'),
('20210226200532'),
('20210226202745'),
('20210226203437'),
('20210226204550'),
('20210226205517'),
('20210301143130'),
('20210301181827');


