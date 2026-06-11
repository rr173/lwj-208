import os
import asyncio
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.database import engine, Base
from app.api import reference, alignment, batch, stats, sample
from app.api.websocket import router as ws_router
from app.sample_data import init_sample_data
from app.services.batch_service import task_manager

Base.metadata.create_all(bind=engine)

app = FastAPI(
    title="Gene Sequence Alignment & Variant Annotation Service",
    description="DNA sequence alignment with Smith-Waterman algorithm, variant annotation, and population analysis",
    version="1.1.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(reference.router)
app.include_router(alignment.router)
app.include_router(batch.router)
app.include_router(stats.router)
app.include_router(sample.router)
app.include_router(ws_router)


@app.on_event("startup")
async def startup_event():
    task_manager.set_loop(asyncio.get_running_loop())
    if os.getenv("INIT_SAMPLE_DATA", "false").lower() == "true":
        init_sample_data()


@app.get("/")
def root():
    return {
        "message": "Gene Sequence Alignment & Variant Annotation Service",
        "version": "1.0.0",
        "docs": "/docs",
    }


@app.get("/health")
def health_check():
    return {"status": "healthy"}
