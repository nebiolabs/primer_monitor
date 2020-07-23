from sqlalchemy.orm import Session

from . import models, schemas


def get_primer(db: Session, primer_id: int):
    return db.query(models.Primer).filter(models.Primer.id == primer_id).first()


def get_primer_by_name(db: Session, name: str):
    return db.query(models.Primer).filter(models.Primer.name == name).first()


def get_primers(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Primer).offset(skip).limit(limit).all()


def create_primer(db: Session, primer: schemas.PrimerCreate):
    db_primer = models.Primer(
        name=primer.name, 
        sequence=primer.sequence, 
        expected_genome=primer.expected_genome,
        primer_set_id=primer.primer_set_id
        )
    db.add(db_primer)
    db.commit()
    db.refresh(db_primer)
    return db_primer


def get_primer_set(db: Session, primer_set_id: int):
    return db.query(models.Primer_set).filter(models.Primer_set.id == primer_set_id).first()


def get_primer_set_by_name(db: Session, name: str):
    return db.query(models.Primer_set).filter(models.Primer_set.name == name).first()


def get_primer_sets(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Primer_set).offset(skip).limit(limit).all()


def create_primer_set(db: Session, primer_set: schemas.PrimerSetCreate):
    db_primer_set = models.Primer_set(
        name=primer_set.name
        )
    db.add(db_primer_set)
    db.commit()
    db.refresh(db_primer_set)
    return db_primer_set


def get_organism(db: Session, organism_id: int):
    return db.query(models.Organism).filter(models.Organism.id == organism_id).first()


def get_organism_by_name(db: Session, name: str):
    return db.query(models.Organism).filter(models.Organism.name == name).first()


def get_organisms(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Organism).offset(skip).limit(limit).all()


def create_organism(db: Session, organism: schemas.OrganismCreate):
    db_organism = models.Organism(
        name=organism.name,
        ncbi_id=organism.ncbi_id,
        always_show=organism.always_show
        )
    db.add(db_organism)
    db.commit()
    db.refresh(db_organism)
    return db_organism


def get_blast_result(db: Session, blast_result_id: int):
    return db.query(models.Blast_result).filter(models.Blast_result.id == blast_result_id).first()


def get_blast_result_by_name(db: Session, name: str):
    return db.query(models.Blast_result).filter(models.Blast_result.name == name).first()


def get_blast_results(db: Session, skip: int = 0, limit: int = 100):
    return db.query(models.Blast_result).offset(skip).limit(limit).all()


def create_blast_result(db: Session, blast_result: schemas.Blast_resultCreate):
    db_blast_result = models.Blast_result(
        name=blast_result.name,
        num_identities=blast_result.num_identities,
        primer_id=blast_result.primer_id,
        organism_id=blast_result.organism_id
        )
    db.add(db_blast_result)
    db.commit()
    db.refresh(db_blast_result)
    return db_blast_result