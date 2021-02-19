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
    geo_location_id integer,
    variant_name character varying
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
-- Name: geo_locations; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.geo_locations (
    id bigint NOT NULL,
    parent_location character varying,
    region character varying NOT NULL,
    division character varying NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
);


--
-- Name: geo_locations_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.geo_locations_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: geo_locations_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.geo_locations_id_seq OWNED BY public.geo_locations.id;


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
    category public.oligo_category
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
            geo_locations.region,
            geo_locations.division,
            geo_locations.id AS geo_id,
            fasta_records.date_collected
           FROM (((public.variant_sites
             JOIN public.fasta_records ON ((variant_sites.fasta_record_id = fasta_records.id)))
             JOIN public.geo_locations ON ((fasta_records.geo_location_id = geo_locations.id)))
             JOIN public.oligos ON ((((variant_sites.ref_start >= oligos.ref_start) AND (variant_sites.ref_start <= oligos.ref_end)) OR ((variant_sites.ref_end >= oligos.ref_start) AND (variant_sites.ref_end <= oligos.ref_end)) OR ((variant_sites.ref_start < oligos.ref_start) AND (variant_sites.ref_end > oligos.ref_end)))))
          WHERE ((((variant_sites.variant_type)::text = 'D'::text) OR ((variant_sites.variant_type)::text = 'X'::text)) AND ((variant_sites.variant)::text !~~ '%N%'::text))
        ), count_query AS (
         SELECT count(*) AS geo_id_count,
            geo_locations.id,
            geo_locations.region,
            geo_locations.division
           FROM (public.fasta_records
             JOIN public.geo_locations ON ((fasta_records.geo_location_id = geo_locations.id)))
          GROUP BY geo_locations.id, geo_locations.region, geo_locations.division
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
            geo_locations.region,
            geo_locations.division,
            geo_locations.id AS geo_id,
            fasta_records.date_collected
           FROM (((public.variant_sites
             JOIN public.fasta_records ON ((variant_sites.fasta_record_id = fasta_records.id)))
             JOIN public.geo_locations ON ((fasta_records.geo_location_id = geo_locations.id)))
             JOIN public.oligos ON ((((variant_sites.ref_start >= oligos.ref_start) AND (variant_sites.ref_start <= oligos.ref_end)) OR ((variant_sites.ref_end >= oligos.ref_start) AND (variant_sites.ref_end <= oligos.ref_end)) OR ((variant_sites.ref_start < oligos.ref_start) AND (variant_sites.ref_end > oligos.ref_end)))))
          WHERE (((variant_sites.variant_type)::text = 'I'::text) AND ((variant_sites.variant)::text !~~ '%N%'::text))
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
    big_query.division,
    big_query.geo_id,
    big_query.date_collected,
    count_query.geo_id_count,
    generate_series(lower((numrange((coord_overlaps.oligo_start)::numeric, (coord_overlaps.oligo_end)::numeric) * numrange((coord_overlaps.variant_start)::numeric, (coord_overlaps.variant_end)::numeric))), (upper((numrange((coord_overlaps.oligo_start)::numeric, (coord_overlaps.oligo_end)::numeric) * numrange((coord_overlaps.variant_start)::numeric, (coord_overlaps.variant_end)::numeric))) - (1)::numeric)) AS coords
   FROM ((big_query coord_overlaps
     JOIN big_query ON (((coord_overlaps.oligo_id = big_query.oligo_id) AND (coord_overlaps.variant_id = big_query.variant_id))))
     JOIN count_query ON ((count_query.id = big_query.geo_id)))
UNION
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
    insert_query.division,
    insert_query.geo_id,
    insert_query.date_collected,
    count_query.geo_id_count,
    insert_query.variant_start AS coords
   FROM ((insert_query coord_overlaps
     JOIN insert_query ON (((coord_overlaps.oligo_id = insert_query.oligo_id) AND (coord_overlaps.variant_id = insert_query.variant_id))))
     JOIN count_query ON ((count_query.id = insert_query.geo_id)))
  WITH NO DATA;


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
    updated_at timestamp(6) without time zone NOT NULL
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
-- Name: primer_sets; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.primer_sets (
    id bigint NOT NULL,
    name character varying,
    user_id bigint,
    organism_id integer NOT NULL,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    status public.primer_set_status DEFAULT 'pending'::public.primer_set_status
);


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
    uid character varying
);


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
-- Name: blast_hits id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.blast_hits ALTER COLUMN id SET DEFAULT nextval('public.blast_hits_id_seq'::regclass);


--
-- Name: fasta_records id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.fasta_records ALTER COLUMN id SET DEFAULT nextval('public.fasta_records_id_seq'::regclass);


--
-- Name: geo_locations id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.geo_locations ALTER COLUMN id SET DEFAULT nextval('public.geo_locations_id_seq'::regclass);


--
-- Name: oligos id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.oligos ALTER COLUMN id SET DEFAULT nextval('public.oligos_id_seq'::regclass);


--
-- Name: organisms id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.organisms ALTER COLUMN id SET DEFAULT nextval('public.organisms_id_seq'::regclass);


--
-- Name: primer_sets id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.primer_sets ALTER COLUMN id SET DEFAULT nextval('public.primer_sets_id_seq'::regclass);


--
-- Name: roles id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.roles ALTER COLUMN id SET DEFAULT nextval('public.roles_id_seq'::regclass);


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
-- Name: fasta_records fasta_records_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.fasta_records
    ADD CONSTRAINT fasta_records_pkey PRIMARY KEY (id);


--
-- Name: geo_locations geo_locations_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.geo_locations
    ADD CONSTRAINT geo_locations_pkey PRIMARY KEY (id);


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
-- Name: primer_sets primer_sets_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.primer_sets
    ADD CONSTRAINT primer_sets_pkey PRIMARY KEY (id);


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
-- Name: index_blast_hits_on_oligo_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_blast_hits_on_oligo_id ON public.blast_hits USING btree (oligo_id);


--
-- Name: index_blast_hits_on_organism_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_blast_hits_on_organism_id ON public.blast_hits USING btree (organism_id);


--
-- Name: index_fasta_records_on_geo_location_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_fasta_records_on_geo_location_id ON public.fasta_records USING btree (geo_location_id);


--
-- Name: index_fasta_records_on_strain; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX index_fasta_records_on_strain ON public.fasta_records USING btree (strain);


--
-- Name: index_oligos_on_primer_set_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_oligos_on_primer_set_id ON public.oligos USING btree (primer_set_id);


--
-- Name: index_primer_sets_on_organism_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_primer_sets_on_organism_id ON public.primer_sets USING btree (organism_id);


--
-- Name: index_primer_sets_on_user_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_primer_sets_on_user_id ON public.primer_sets USING btree (user_id);


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
-- Name: user_roles fk_rails_3369e0d5fc; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.user_roles
    ADD CONSTRAINT fk_rails_3369e0d5fc FOREIGN KEY (role_id) REFERENCES public.roles(id);


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
-- Name: fasta_records fk_rails_8783be51d9; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.fasta_records
    ADD CONSTRAINT fk_rails_8783be51d9 FOREIGN KEY (geo_location_id) REFERENCES public.geo_locations(id);


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
('20210219153745'),
('20210219185233');


