from typing import List

from fastapi import Depends, APIRouter, HTTPException
from sqlalchemy.orm import Session

from primer_monitor_api import crud, schemas
from primer_monitor_api.api.deps import get_db

router = APIRouter()

@router.post("/", response_model=schemas.PrimerSet)
def create_primer_set(primer_set: schemas.PrimerSetCreate, db: Session = Depends(get_db)):
    db_primer_set = crud.get_primer_set_by_name(db, name=primer_set.name)
    if db_primer_set:
        raise HTTPException(status_code=400, detail="PrimerSet name already registered")
    return crud.create_primer_set(db=db, primer_set=primer_set)


@router.get("/", response_model=List[schemas.PrimerSet])
def read_primer_sets(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    primer_sets = crud.get_primer_sets(db, skip=skip, limit=limit)
    return primer_sets


@router.get("/{primer_set_id}", response_model=schemas.PrimerSet)
def read_primer_set(primer_set_id: int, db: Session = Depends(get_db)):
    db_primer_set = crud.get_primer_set(db, primer_set_id=primer_set_id)
    if db_primer_set is None:
        raise HTTPException(status_code=404, detail="Primer set not found")
    return db_primer_set

