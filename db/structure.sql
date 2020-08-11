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
-- Name: amplicons; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.amplicons (
    id bigint NOT NULL,
    name character varying,
    user_id bigint,
    organism_id bigint,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
);


--
-- Name: amplicons_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE public.amplicons_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: amplicons_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE public.amplicons_id_seq OWNED BY public.amplicons.id;


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
-- Name: oligos; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.oligos (
    id bigint NOT NULL,
    name character varying,
    sequence character varying,
    amplicon_id bigint,
    created_at timestamp(6) without time zone NOT NULL,
    updated_at timestamp(6) without time zone NOT NULL
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
-- Name: schema_migrations; Type: TABLE; Schema: public; Owner: -
--

CREATE TABLE public.schema_migrations (
    version character varying NOT NULL
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
    updated_at timestamp(6) without time zone NOT NULL
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
-- Name: amplicons id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.amplicons ALTER COLUMN id SET DEFAULT nextval('public.amplicons_id_seq'::regclass);


--
-- Name: blast_hits id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.blast_hits ALTER COLUMN id SET DEFAULT nextval('public.blast_hits_id_seq'::regclass);


--
-- Name: oligos id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.oligos ALTER COLUMN id SET DEFAULT nextval('public.oligos_id_seq'::regclass);


--
-- Name: organisms id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.organisms ALTER COLUMN id SET DEFAULT nextval('public.organisms_id_seq'::regclass);


--
-- Name: users id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.users ALTER COLUMN id SET DEFAULT nextval('public.users_id_seq'::regclass);


--
-- Name: amplicons amplicons_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.amplicons
    ADD CONSTRAINT amplicons_pkey PRIMARY KEY (id);


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
-- Name: schema_migrations schema_migrations_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.schema_migrations
    ADD CONSTRAINT schema_migrations_pkey PRIMARY KEY (version);


--
-- Name: users users_pkey; Type: CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.users
    ADD CONSTRAINT users_pkey PRIMARY KEY (id);


--
-- Name: index_amplicons_on_organism_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_amplicons_on_organism_id ON public.amplicons USING btree (organism_id);


--
-- Name: index_amplicons_on_user_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_amplicons_on_user_id ON public.amplicons USING btree (user_id);


--
-- Name: index_blast_hits_on_oligo_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_blast_hits_on_oligo_id ON public.blast_hits USING btree (oligo_id);


--
-- Name: index_blast_hits_on_organism_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_blast_hits_on_organism_id ON public.blast_hits USING btree (organism_id);


--
-- Name: index_oligos_on_amplicon_id; Type: INDEX; Schema: public; Owner: -
--

CREATE INDEX index_oligos_on_amplicon_id ON public.oligos USING btree (amplicon_id);


--
-- Name: blast_hits fk_rails_1f04a34db0; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.blast_hits
    ADD CONSTRAINT fk_rails_1f04a34db0 FOREIGN KEY (organism_id) REFERENCES public.organisms(id);


--
-- Name: amplicons fk_rails_52b9cf4012; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.amplicons
    ADD CONSTRAINT fk_rails_52b9cf4012 FOREIGN KEY (organism_id) REFERENCES public.organisms(id);


--
-- Name: oligos fk_rails_58e7536518; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.oligos
    ADD CONSTRAINT fk_rails_58e7536518 FOREIGN KEY (amplicon_id) REFERENCES public.amplicons(id);


--
-- Name: amplicons fk_rails_a78cff2c70; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY public.amplicons
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
('20200810101557');


