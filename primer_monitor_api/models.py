from sqlalchemy import Boolean, Column, ForeignKey, Integer, String
from sqlalchemy.orm import relationship

from .database import Base


class Primer_set(Base):
    # CREATE TABLE public.primer_sets (
    # 	id serial UNIQUE NOT NULL,
    # 	"name" varchar UNIQUE NOT NULL
    # );
    __tablename__ = "primer_sets"

    id = Column(Integer, primary_key=True, index=True,)
    name = Column(String, unique=True, index=True,)
    primers = relationship("Primer", back_populates="primer_set")


class Primer(Base):
    # CREATE TABLE public.primers (
    # 	id serial NOT NULL,
    # 	"name" varchar NOT NULL,
    # 	"sequence" varchar NOT NULL,
    #         expected_genome varchar NOT NULL,
    #         primer_set_id int NOT NULL references primer_sets(id),
    # 	CONSTRAINT primers_id_key UNIQUE (id),
    # 	CONSTRAINT primers_name_key UNIQUE ("name"),
    # 	CONSTRAINT primers_sequence_key UNIQUE ("sequence")
    # );
    __tablename__ = "primers"

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String, unique=True, index=True)
    sequence = Column(String, unique=True)
    expected_genome = Column(String, unique=True)

    primer_set_id = Column(Integer, ForeignKey("primer_sets.id"))
    primer_set = relationship("Primer_set", back_populates="primers")
    
    blast_results = relationship("Blast_result", back_populates='primer')


class Organism(Base):
    # CREATE TABLE public.organisms (
    # 	id serial NOT NULL,
    # 	"name" varchar NOT NULL,
    # 	ncbi_id varchar NOT NULL,
    # 	always_show bool NOT NULL DEFAULT false,
    # 	CONSTRAINT organisms_id_key UNIQUE (id),
    # 	CONSTRAINT organisms_name_key UNIQUE ("name"),
    # 	CONSTRAINT organisms_ncbi_id_key UNIQUE (ncbi_id)
    # );
    __tablename__ = "organisms"

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String, unique=True, index=True)
    ncbi_id = Column(String, unique=True)
    always_show = Column(Boolean, default=False)
    blast_results = relationship("Blast_result", back_populates='organism')



class Blast_result(Base):
    # CREATE TABLE public.blast_results (
    # 	id serial unique NOT null primary key,
    # 	primer_id int NOT null references primers(id), 
    # 	organism_id int NOT null references organisms(id),
    # 	num_identities int not null
        
    # );
    __tablename__ = "blast_results"

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String, unique=True, index=True)
    num_identities = Column(Integer)

    primer_id = Column(Integer, ForeignKey("primers.id"))
    primer = relationship("Primer", back_populates='blast_results')

    organism_id = Column(Integer, ForeignKey("organisms.id"))
    organism = relationship("Organism", back_populates='blast_results')