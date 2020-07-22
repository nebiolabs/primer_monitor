
CREATE TABLE public.primer_sets (
	id serial UNIQUE NOT NULL,
	"name" varchar UNIQUE NOT NULL
);

CREATE TABLE public.primers (
	id serial NOT NULL,
	"name" varchar NOT NULL,
	"sequence" varchar NOT NULL,
        expected_genome varchar NOT NULL,
        primer_set_id int NOT NULL references primer_sets(id),
	CONSTRAINT primers_id_key UNIQUE (id),
	CONSTRAINT primers_name_key UNIQUE ("name"),
	CONSTRAINT primers_sequence_key UNIQUE ("sequence")
);

CREATE TABLE public.organisms (
	id serial NOT NULL,
	"name" varchar NOT NULL,
	ncbi_id varchar NOT NULL,
	always_show bool NOT NULL DEFAULT false,
	CONSTRAINT organisms_id_key UNIQUE (id),
	CONSTRAINT organisms_name_key UNIQUE ("name"),
	CONSTRAINT organisms_ncbi_id_key UNIQUE (ncbi_id)
);

CREATE TABLE public.blast_results (
	id serial unique NOT null primary key,
	primer_id int NOT null references primers(id), 
	organism_id int NOT null references organisms(id),
	num_identities int not null
	
);

