from typing import List, Optional

from pydantic import BaseModel


from sqlalchemy import Boolean, Column, ForeignKey, Integer, String
from sqlalchemy.orm import relationship

from .database import Base



class PrimerBase(BaseModel):
    name: str
    sequence: str
    expected_genome: str
    primer_set_id: int


class PrimerCreate(PrimerBase):
    pass


class Primer(PrimerBase):
    id: int 

    class Config:
        orm_mode = True


class PrimerSetBase(BaseModel):
    name: str


class PrimerSetCreate(PrimerSetBase):
    # ie passwords
    pass


class PrimerSet(PrimerSetBase):
    id: int 
    primers: List[Primer] = []

    class Config:
        orm_mode = True


class OrganismBase(BaseModel):
    name: str
    ncbi_id: str
    always_show: bool


class OrganismCreate(OrganismBase):
    pass


class Organism(OrganismBase):
    id: int 

    class Config:
        orm_mode = True


class Blast_resultBase(BaseModel):
    name: str
    num_identities: int

# Properties to receive via API on creation
class Blast_resultCreate(Blast_resultBase):
    primer_id: int
    organism_id: int

class Blast_result(Blast_resultBase):
    id: int 
    primer: Primer
    organism: Organism

    class Config:
        orm_mode = True