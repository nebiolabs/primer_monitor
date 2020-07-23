from typing import List

from fastapi import Depends, APIRouter, HTTPException
from sqlalchemy.orm import Session

from primer_monitor_api import crud, schemas
from primer_monitor_api.api.deps import get_db

router = APIRouter()

@router.post("/", response_model=schemas.Primer)
def create_primer(primer: schemas.PrimerCreate, db: Session = Depends(get_db)):
    db_primer = crud.get_primer_by_name(db, name=primer.name)
    if db_primer:
        raise HTTPException(status_code=400, detail="Primer name already registered")
    return crud.create_primer(db=db, primer=primer)


@router.get("/", response_model=List[schemas.Primer])
def read_primers(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    primers = crud.get_primers(db, skip=skip, limit=limit)
    return primers


@router.get("/{primer_id}", response_model=schemas.Primer)
def read_primer(primer_id: int, db: Session = Depends(get_db)):
    db_primer = crud.get_primer(db, primer_id=primer_id)
    if db_primer is None:
        raise HTTPException(status_code=404, detail="Primer not found")
    return db_primer
