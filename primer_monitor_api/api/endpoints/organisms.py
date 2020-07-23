from typing import List

from fastapi import Depends, APIRouter, HTTPException
from sqlalchemy.orm import Session

from primer_monitor_api import crud, schemas
from primer_monitor_api.api.deps import get_db

router = APIRouter()

@router.post("/", response_model=schemas.Organism)
def create_organism(organism: schemas.OrganismCreate, db: Session = Depends(get_db)):
    db_organism = crud.get_organism_by_name(db, name=organism.name)
    if db_organism:
        raise HTTPException(status_code=400, detail="Organism name already registered")
    return crud.create_organism(db=db, organism=organism)


@router.get("/", response_model=List[schemas.Organism])
def read_organisms(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    organisms = crud.get_organisms(db, skip=skip, limit=limit)
    return organisms


@router.get("/{organism_id}", response_model=schemas.Organism)
def read_organism(organism_id: int, db: Session = Depends(get_db)):
    db_organism = crud.get_organism(db, organism_id=organism_id)
    if db_organism is None:
        raise HTTPException(status_code=404, detail="Primer set not found")
    return db_organism

