from typing import List

from fastapi import Depends, APIRouter, HTTPException
from sqlalchemy.orm import Session

from primer_monitor_api import crud, schemas
from primer_monitor_api.api.deps import get_db

router = APIRouter()

@router.post("/", response_model=schemas.Blast_result)
def create_blast_result(blast_result: schemas.Blast_resultCreate, db: Session = Depends(get_db)):
    db_blast_result = crud.get_blast_result_by_name(db, name=blast_result.name)
    if db_blast_result:
        raise HTTPException(status_code=400, detail="Blast_result name already registered")
    return crud.create_blast_result(db=db, blast_result=blast_result)


@router.get("/", response_model=List[schemas.Blast_result])
def read_blast_results(skip: int = 0, limit: int = 100, db: Session = Depends(get_db)):
    blast_results = crud.get_blast_results(db, skip=skip, limit=limit)
    return blast_results


@router.get("/{blast_result_id}", response_model=schemas.Blast_result)
def read_blast_result(blast_result_id: int, db: Session = Depends(get_db)):
    db_blast_result = crud.get_blast_result(db, blast_result_id=blast_result_id)
    if db_blast_result is None:
        raise HTTPException(status_code=404, detail="Primer set not found")
    return db_blast_result

