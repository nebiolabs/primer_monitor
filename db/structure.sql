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
    updated_at timestamp(6) without time zone NOT NULL,
    detailed_geo_location_alias_id bigint NOT NULL
);


--
-- Name: fasta_records; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.fasta_records (
    id bigint NOT NULL,
    strain character varying NOT NULL,
    genbank_accession character varying,
    gisaid_epi_isl character varying,
    date_collected date,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    variant_name character varying,
    detailed_geo_location_id bigint NOT NULL
);


--
-- Name: counts; Type: MATERIALIZED VIEW; Schema: public; Owner: -
--

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
  WITH NO DATA;


--
-- Name: date_counts; Type: VIEW; Schema: public; Owner: -
--

CREATE VIEW public.date_counts AS
 SELECT count(*) AS total_count,
    fr.detailed_geo_location_id,
    fr.date_collected
   FROM public.fasta_records fr
  GROUP BY fr.detailed_geo_location_id, fr.date_collected;


--
-- Name: detailed_geo_location_aliases; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.detailed_geo_location_aliases (
    id bigint NOT NULL,
    world character varying NOT NULL,
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
    detailed_geo_location_alias_id bigint NOT NULL
);


--
-- Name: join_subscribed_location_to_ids; Type: VIEW; Schema: public; Owner: -
--

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


--
-- Name: oligos; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.oligos (
    id bigint NOT NULL,
    name character varying,
    sequence character varying,
    primer_set_id bigint,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    locus character varying,
    category public.oligo_category,
    ref_start bigint,
    ref_end bigint,
    short_name character varying,
    strand character varying
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
-- Name: time_counts; Type: MATERIALIZED VIEW; Schema: public; Owner: -
--

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
  WITH NO DATA;


--
-- Name: variant_sites; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.variant_sites (
    id bigint NOT NULL,
    ref_start integer NOT NULL,
    variant_type character varying NOT NULL,
    variant character varying NOT NULL,
    fasta_record_id bigint NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    ref_end integer NOT NULL
);


--
-- Name: variant_overlaps; Type: MATERIALIZED VIEW; Schema: public; Owner: -
--

CREATE MATERIALIZED VIEW public.variant_overlaps AS
 SELECT oligos.id AS oligo_id,
    oligos.name AS oligo_name,
    oligos.ref_start AS oligo_start,
    oligos.ref_end AS oligo_end,
    oligos.short_name AS oligo_short_name,
    primer_sets.id AS primer_set_id,
    primer_sets.name AS primer_set_name,
    variant_sites.id AS variant_id,
    variant_sites.variant_type,
    variant_sites.variant,
    variant_sites.ref_start AS variant_start,
    variant_sites.ref_end AS variant_end,
    COALESCE(detailed_geo_locations.region, ''::character varying) AS region,
    COALESCE(detailed_geo_locations.subregion, ''::character varying) AS subregion,
    COALESCE(detailed_geo_locations.division, ''::character varying) AS division,
    COALESCE(detailed_geo_locations.subdivision, ''::character varying) AS subdivision,
    detailed_geo_locations.id AS detailed_geo_location_id,
    COALESCE(fasta_records.date_collected, '1900-01-01'::date) AS date_collected
   FROM ((((public.variant_sites
     JOIN public.fasta_records ON ((variant_sites.fasta_record_id = fasta_records.id)))
     JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
     JOIN public.oligos ON ((((variant_sites.ref_start >= oligos.ref_start) AND (variant_sites.ref_start < oligos.ref_end)) OR ((variant_sites.ref_end > oligos.ref_start) AND (variant_sites.ref_end <= oligos.ref_end)) OR ((variant_sites.ref_start < oligos.ref_start) AND (variant_sites.ref_end > oligos.ref_end)))))
     JOIN public.primer_sets ON ((oligos.primer_set_id = primer_sets.id)))
  WHERE ((((variant_sites.variant_type)::text = 'D'::text) OR ((variant_sites.variant_type)::text = 'X'::text)) AND ((variant_sites.variant)::text !~~ '%N%'::text))
UNION ALL
 SELECT oligos.id AS oligo_id,
    oligos.name AS oligo_name,
    oligos.ref_start AS oligo_start,
    oligos.ref_end AS oligo_end,
    oligos.short_name AS oligo_short_name,
    primer_sets.id AS primer_set_id,
    primer_sets.name AS primer_set_name,
    variant_sites.id AS variant_id,
    variant_sites.variant_type,
    variant_sites.variant,
    variant_sites.ref_start AS variant_start,
    variant_sites.ref_end AS variant_end,
    COALESCE(detailed_geo_locations.region, ''::character varying) AS region,
    COALESCE(detailed_geo_locations.subregion, ''::character varying) AS subregion,
    COALESCE(detailed_geo_locations.division, ''::character varying) AS division,
    COALESCE(detailed_geo_locations.subdivision, ''::character varying) AS subdivision,
    detailed_geo_locations.id AS detailed_geo_location_id,
    COALESCE(fasta_records.date_collected, '1900-01-01'::date) AS date_collected
   FROM ((((public.variant_sites
     JOIN public.fasta_records ON ((variant_sites.fasta_record_id = fasta_records.id)))
     JOIN public.detailed_geo_locations ON ((fasta_records.detailed_geo_location_id = detailed_geo_locations.id)))
     JOIN public.oligos ON (((variant_sites.ref_start > oligos.ref_start) AND (variant_sites.ref_start <= oligos.ref_end))))
     JOIN public.primer_sets ON ((oligos.primer_set_id = primer_sets.id)))
  WHERE (((variant_sites.variant_type)::text = 'I'::text) AND ((variant_sites.variant)::text !~~ '%N%'::text))
  WITH NO DATA;


--
-- Name: oligo_variant_overlaps; Type: MATERIALIZED VIEW; Schema: public; Owner: -
--

CREATE MATERIALIZED VIEW public.oligo_variant_overlaps AS
 SELECT variant_overlaps.oligo_id,
    variant_overlaps.oligo_name,
    variant_overlaps.oligo_start,
    variant_overlaps.oligo_end,
    variant_overlaps.oligo_short_name,
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
    counts.region_count,
    counts.region_subregion_count,
    counts.region_subregion_division_count,
    counts.region_subregion_division_subdivision_count,
    time_counts.region_time_count,
    time_counts.region_subregion_time_count,
    time_counts.region_subregion_division_time_count,
    time_counts.region_subregion_division_subdivision_time_count,
    generate_series(lower((numrange((variant_overlaps.oligo_start)::numeric, (variant_overlaps.oligo_end)::numeric) * numrange((variant_overlaps.variant_start)::numeric, (variant_overlaps.variant_end)::numeric))), (upper((numrange((variant_overlaps.oligo_start)::numeric, (variant_overlaps.oligo_end)::numeric) * numrange((variant_overlaps.variant_start)::numeric, (variant_overlaps.variant_end)::numeric))) - (1)::numeric)) AS coords
   FROM ((public.variant_overlaps variant_overlaps
     JOIN public.counts ON ((((counts.region)::text = (variant_overlaps.region)::text) AND ((counts.subregion)::text = (variant_overlaps.subregion)::text) AND ((counts.division)::text = (variant_overlaps.division)::text) AND ((counts.subdivision)::text = (variant_overlaps.subdivision)::text))))
     JOIN public.time_counts ON ((((time_counts.region)::text = (variant_overlaps.region)::text) AND ((time_counts.subregion)::text = (variant_overlaps.subregion)::text) AND ((time_counts.division)::text = (variant_overlaps.division)::text) AND ((time_counts.subdivision)::text = (variant_overlaps.subdivision)::text) AND (time_counts.date_collected = variant_overlaps.date_collected))))
  WHERE ((variant_overlaps.variant_type)::text <> 'I'::text)
UNION ALL
 SELECT variant_overlaps.oligo_id,
    variant_overlaps.oligo_name,
    variant_overlaps.oligo_start,
    variant_overlaps.oligo_end,
    variant_overlaps.oligo_short_name,
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
    counts.region_count,
    counts.region_subregion_count,
    counts.region_subregion_division_count,
    counts.region_subregion_division_subdivision_count,
    time_counts.region_time_count,
    time_counts.region_subregion_time_count,
    time_counts.region_subregion_division_time_count,
    time_counts.region_subregion_division_subdivision_time_count,
    variant_overlaps.variant_start AS coords
   FROM ((public.variant_overlaps
     JOIN public.counts ON ((((counts.region)::text = (variant_overlaps.region)::text) AND ((counts.subregion)::text = (variant_overlaps.subregion)::text) AND ((counts.division)::text = (variant_overlaps.division)::text) AND ((counts.subdivision)::text = (variant_overlaps.subdivision)::text))))
     JOIN public.time_counts ON ((((time_counts.region)::text = (variant_overlaps.region)::text) AND ((time_counts.subregion)::text = (variant_overlaps.subregion)::text) AND ((time_counts.division)::text = (variant_overlaps.division)::text) AND ((time_counts.subdivision)::text = (variant_overlaps.subdivision)::text) AND (time_counts.date_collected = variant_overlaps.date_collected))))
  WHERE ((variant_overlaps.variant_type)::text = 'I'::text)
  WITH NO DATA;


--
-- Name: primer_set_subscriptions; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.primer_set_subscriptions (
    id bigint NOT NULL,
    primer_set_id bigint NOT NULL,
    user_id bigint NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
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
-- Name: identify_primers_for_notifications; Type: MATERIALIZED VIEW; Schema: public; Owner: -
--

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
  WITH NO DATA;


--
-- Name: initial_score; Type: MATERIALIZED VIEW; Schema: public; Owner: -
--

CREATE MATERIALIZED VIEW public.initial_score AS
 WITH all_combos AS (
         SELECT oligos_1.oligo_name,
            oligos_1.locus,
            oligos_1.primer_set_id,
            fa_1.detailed_geo_location_id,
            fa_1.date_collected
           FROM (( SELECT DISTINCT oligos.name AS oligo_name,
                    oligos.locus,
                    oligos.primer_set_id
                   FROM public.oligos) oligos_1
             CROSS JOIN ( SELECT DISTINCT fasta_records.detailed_geo_location_id,
                    fasta_records.date_collected
                   FROM public.fasta_records
                  WHERE (fasta_records.date_collected >= '2020-01-01'::date)) fa_1)
        ), variants AS (
         SELECT count(variant_sites.fasta_record_id) AS observed_count,
            ovo.detailed_geo_location_id,
            ((count(*))::double precision / ((5)::double precision * (count(variant_sites.fasta_record_id))::double precision)) AS three_p_score,
            oligos.locus,
            ovo.primer_set_id,
            ovo.oligo_name,
            ovo.date_collected
           FROM (((public.oligo_variant_overlaps ovo
             JOIN public.variant_sites ON ((variant_sites.id = ovo.variant_id)))
             JOIN public.fasta_records ON ((fasta_records.id = variant_sites.fasta_record_id)))
             JOIN public.oligos ON (((oligos.primer_set_id = ovo.primer_set_id) AND ((oligos.name)::text = (ovo.oligo_name)::text))))
          WHERE (((ovo.oligo_end)::numeric - ovo.coords) <= (5)::numeric)
          GROUP BY ovo.detailed_geo_location_id, oligos.locus, ovo.primer_set_id, ovo.oligo_name, ovo.date_collected
        )
 SELECT all_combos.date_collected,
    all_combos.detailed_geo_location_id,
    COALESCE(variants.observed_count, (0)::bigint) AS observed_count,
    COALESCE(variants.three_p_score, (0)::double precision) AS three_p_score,
    all_combos.locus,
    all_combos.primer_set_id,
    all_combos.oligo_name
   FROM (all_combos
     LEFT JOIN variants ON (((all_combos.date_collected = variants.date_collected) AND (all_combos.primer_set_id = variants.primer_set_id) AND ((all_combos.oligo_name)::text = (variants.oligo_name)::text) AND (all_combos.detailed_geo_location_id = variants.detailed_geo_location_id))))
  WITH NO DATA;


--
-- Name: join_subscribed_location_to_id; Type: VIEW; Schema: public; Owner: -
--

CREATE VIEW public.join_subscribed_location_to_id AS
 WITH subscribed_ids AS (
         SELECT subscribed_geo_locations.user_id,
            subscribed_geo_locations.detailed_geo_location_alias_id,
            detailed_geo_location_aliases_1.region,
            detailed_geo_location_aliases_1.subregion,
            detailed_geo_location_aliases_1.division,
            detailed_geo_location_aliases_1.subdivision
           FROM (public.detailed_geo_location_aliases detailed_geo_location_aliases_1
             JOIN public.subscribed_geo_locations ON ((subscribed_geo_locations.detailed_geo_location_alias_id = detailed_geo_location_aliases_1.id)))
        )
 SELECT subscribed_ids.user_id,
    subscribed_ids.detailed_geo_location_alias_id AS subscribed_id,
    detailed_geo_location_aliases.id AS detailed_geo_location_alias_id,
    detailed_geo_locations.id AS detailed_geo_location_id
   FROM ((public.detailed_geo_location_aliases
     JOIN subscribed_ids ON ((((subscribed_ids.region IS NULL) OR ((subscribed_ids.region)::text = (detailed_geo_location_aliases.region)::text)) AND ((subscribed_ids.subregion IS NULL) OR ((subscribed_ids.subregion)::text = (detailed_geo_location_aliases.subregion)::text)) AND ((subscribed_ids.division IS NULL) OR ((subscribed_ids.division)::text = (detailed_geo_location_aliases.division)::text)) AND ((subscribed_ids.subdivision IS NULL) OR ((subscribed_ids.subdivision)::text = (detailed_geo_location_aliases.subdivision)::text)))))
     JOIN public.detailed_geo_locations ON ((detailed_geo_locations.detailed_geo_location_alias_id = detailed_geo_location_aliases.id)));


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
    primer_set_id bigint NOT NULL,
    oligo_id bigint NOT NULL,
    verified_notification_id bigint,
    coordinate integer NOT NULL,
    fraction_variant double precision NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    subscribed_geo_location_id bigint NOT NULL,
    primer_set_subscription_id bigint NOT NULL,
    user_id bigint NOT NULL
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
    user_id bigint NOT NULL,
    date_sent date,
    status character varying NOT NULL,
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
-- Name: counts_region_subregion_division_subdivision_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX counts_region_subregion_division_subdivision_idx ON public.counts USING btree (region, subregion, division, subdivision);


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
-- Name: index_detailed_geo_locations_on_detailed_geo_location_alias_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_detailed_geo_locations_on_detailed_geo_location_alias_id ON public.detailed_geo_locations USING btree (detailed_geo_location_alias_id);


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
-- Name: index_proposed_notifications_on_oligo_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_proposed_notifications_on_oligo_id ON public.proposed_notifications USING btree (oligo_id);


--
-- Name: index_proposed_notifications_on_primer_set_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_proposed_notifications_on_primer_set_id ON public.proposed_notifications USING btree (primer_set_id);


--
-- Name: index_proposed_notifications_on_primer_set_subscription_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_proposed_notifications_on_primer_set_subscription_id ON public.proposed_notifications USING btree (primer_set_subscription_id);


--
-- Name: index_proposed_notifications_on_subscribed_geo_location_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_proposed_notifications_on_subscribed_geo_location_id ON public.proposed_notifications USING btree (subscribed_geo_location_id);


--
-- Name: index_proposed_notifications_on_user_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_proposed_notifications_on_user_id ON public.proposed_notifications USING btree (user_id);


--
-- Name: index_proposed_notifications_on_verified_notification_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_proposed_notifications_on_verified_notification_id ON public.proposed_notifications USING btree (verified_notification_id);


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
-- Name: index_verified_notifications_on_user_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_verified_notifications_on_user_id ON public.verified_notifications USING btree (user_id);


--
-- Name: time_counts_region_subregion_division_subdivision_date_coll_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX time_counts_region_subregion_division_subdivision_date_coll_idx ON public.time_counts USING btree (region, subregion, division, subdivision, date_collected);


--
-- Name: tmp; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX tmp ON public.subscribed_geo_locations USING btree (detailed_geo_location_alias_id);


--
-- Name: variant_overlaps_region_subregion_division_subdivision_date_idx; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX variant_overlaps_region_subregion_division_subdivision_date_idx ON public.variant_overlaps USING btree (region, subregion, division, subdivision, date_collected);


--
-- Name: proposed_notifications fk_rails_03fe9a3c07; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.proposed_notifications
    ADD CONSTRAINT fk_rails_03fe9a3c07 FOREIGN KEY (subscribed_geo_location_id) REFERENCES public.subscribed_geo_locations(id);


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
    ADD CONSTRAINT fk_rails_32a2a94275 FOREIGN KEY (primer_set_id) REFERENCES public.primer_sets(id);


--
-- Name: user_roles fk_rails_3369e0d5fc; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.user_roles
    ADD CONSTRAINT fk_rails_3369e0d5fc FOREIGN KEY (role_id) REFERENCES public.roles(id);


--
-- Name: proposed_notifications fk_rails_39f3a62e67; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.proposed_notifications
    ADD CONSTRAINT fk_rails_39f3a62e67 FOREIGN KEY (oligo_id) REFERENCES public.oligos(id);


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
-- Name: genomic_features fk_rails_65e85371d4; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.genomic_features
    ADD CONSTRAINT fk_rails_65e85371d4 FOREIGN KEY (organism_id) REFERENCES public.organisms(id);


--
-- Name: verified_notifications fk_rails_660c55653c; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.verified_notifications
    ADD CONSTRAINT fk_rails_660c55653c FOREIGN KEY (user_id) REFERENCES public.users(id);


--
-- Name: proposed_notifications fk_rails_736e9d2781; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.proposed_notifications
    ADD CONSTRAINT fk_rails_736e9d2781 FOREIGN KEY (verified_notification_id) REFERENCES public.verified_notifications(id);


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
-- Name: detailed_geo_locations fk_rails_8a5906d321; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.detailed_geo_locations
    ADD CONSTRAINT fk_rails_8a5906d321 FOREIGN KEY (detailed_geo_location_alias_id) REFERENCES public.detailed_geo_location_aliases(id);


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
-- Name: proposed_notifications fk_rails_bfe6da46a3; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.proposed_notifications
    ADD CONSTRAINT fk_rails_bfe6da46a3 FOREIGN KEY (user_id) REFERENCES public.users(id);


--
-- Name: primer_set_subscriptions fk_rails_e7701775d5; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.primer_set_subscriptions
    ADD CONSTRAINT fk_rails_e7701775d5 FOREIGN KEY (user_id) REFERENCES public.users(id);


--
-- Name: proposed_notifications fk_rails_f9871ce32b; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.proposed_notifications
    ADD CONSTRAINT fk_rails_f9871ce32b FOREIGN KEY (primer_set_subscription_id) REFERENCES public.primer_set_subscriptions(id);


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
('20210301181827'),
('20210302184915'),
('20210302193440'),
('20210302193832'),
('20210303180741'),
('20210303211641'),
('20210304200709'),
('20210307125905'),
('20210312143921'),
('20210312185232'),
('20210313180917'),
('20210318204840'),
('20210321204308'),
('20210323205921'),
('20210325203319'),
('20210405161255'),
('20210413164024'),
('20210413164039'),
('20210619121250'),
('20210619220014'),
('20210621161928');


