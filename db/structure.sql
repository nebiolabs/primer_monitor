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
    region character varying,
    country character varying,
    division character varying,
    date_submitted date,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
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
-- Name: oligos; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.oligos (
    id bigint NOT NULL,
    name character varying,
    sequence character varying,
    primer_set_id bigint,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL,
    start_pos integer,
    end_pos integer
);


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
    updated_at timestamp(6) without time zone NOT NULL
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
    crypted_password character varying,
    password_salt character varying,
    persistence_token character varying,
    single_access_token character varying,
    perishable_token character varying,
    login_count integer DEFAULT 0 NOT NULL,
    failed_login_count integer DEFAULT 0 NOT NULL,
    last_request_at timestamp without time zone,
    current_login_at timestamp without time zone,
    last_login_at timestamp without time zone,
    current_login_ip character varying,
    last_login_ip character varying,
    active boolean DEFAULT false,
    approved boolean DEFAULT false,
    confirmed boolean DEFAULT false,
    send_primer_updates boolean DEFAULT false NOT NULL
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
-- Name: variant_sites; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.variant_sites (
    id bigint NOT NULL,
    "position" integer,
    variant_type character varying,
    variant character varying,
    fasta_record_id bigint,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
);


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
-- Name: index_users_on_perishable_token; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX index_users_on_perishable_token ON public.users USING btree (perishable_token);


--
-- Name: index_users_on_persistence_token; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX index_users_on_persistence_token ON public.users USING btree (persistence_token);


--
-- Name: index_users_on_single_access_token; Type: INDEX; Schema: public; Owner: -
--

CREATE UNIQUE INDEX index_users_on_single_access_token ON public.users USING btree (single_access_token);


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
('20210129031328');


