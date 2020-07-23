from fastapi import APIRouter

from primer_monitor_api.api.endpoints import primer_sets, primers, organisms, blast_results

api_router = APIRouter()
api_router.include_router(primer_sets.router, prefix="/primer_sets", tags=["primer_sets"])
api_router.include_router(primers.router, prefix="/primers", tags=["primers"])
api_router.include_router(organisms.router, prefix="/organisms", tags=["organisms"])
api_router.include_router(blast_results.router, prefix="/blast_results", tags=["blast_results"])