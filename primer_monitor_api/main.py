from typing import List

from fastapi import Depends, FastAPI, HTTPException
from sqlalchemy.orm import Session

from primer_monitor_api import models
from primer_monitor_api.api.api import api_router

from primer_monitor_api.database import engine

models.Base.metadata.create_all(bind=engine)

app = FastAPI(title='Primer Monitor')

app.include_router(api_router)

