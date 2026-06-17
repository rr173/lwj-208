import os
import asyncio
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.database import engine, Base, SessionLocal
from app.api import reference, alignment, batch, stats, sample, phylogeny, scoring, ld, primer, transmission, domain, synteny, audit, mutational_signature, pipeline
from app.api.websocket import router as ws_router
from app.sample_data import init_sample_data
from app.services.batch_service import task_manager
from app.services.scoring_service import seed_default_rules
from app.services.mutational_signature_service import list_reference_signatures
from app.audit import AuditLogMiddleware

Base.metadata.create_all(bind=engine)


def _migrate_pipeline_tables():
    from sqlalchemy import text, inspect
    inspector = inspect(engine)
    with engine.connect() as conn:
        exec_cols = [c["name"] for c in inspector.get_columns("pipeline_executions")]
        if "template_snapshot" not in exec_cols:
            conn.execute(text("ALTER TABLE pipeline_executions ADD COLUMN template_snapshot JSON"))
        if "resume_count" not in exec_cols:
            conn.execute(text("ALTER TABLE pipeline_executions ADD COLUMN resume_count INTEGER DEFAULT 0"))

        step_cols = [c["name"] for c in inspector.get_columns("pipeline_step_executions")]
        if "retry_count" not in step_cols:
            conn.execute(text("ALTER TABLE pipeline_step_executions ADD COLUMN retry_count INTEGER DEFAULT 0"))
        if "last_error" not in step_cols:
            conn.execute(text("ALTER TABLE pipeline_step_executions ADD COLUMN last_error TEXT"))
        if "output_data" not in step_cols:
            conn.execute(text("ALTER TABLE pipeline_step_executions ADD COLUMN output_data JSON"))
        conn.commit()


try:
    _migrate_pipeline_tables()
except Exception:
    pass

app = FastAPI(
    title="Gene Sequence Alignment & Variant Annotation Service",
    description="DNA sequence alignment with Smith-Waterman algorithm, variant annotation, and population analysis",
    version="1.2.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.add_middleware(AuditLogMiddleware)

app.include_router(reference.router)
app.include_router(alignment.router)
app.include_router(batch.router)
app.include_router(stats.router)
app.include_router(sample.router)
app.include_router(phylogeny.router)
app.include_router(scoring.router)
app.include_router(ld.router)
app.include_router(primer.router)
app.include_router(transmission.router)
app.include_router(domain.router)
app.include_router(synteny.router)
app.include_router(audit.router)
app.include_router(ws_router)
app.include_router(mutational_signature.router)
app.include_router(pipeline.router)


@app.on_event("startup")
async def startup_event():
    task_manager.set_loop(asyncio.get_running_loop())
    if os.getenv("INIT_SAMPLE_DATA", "false").lower() == "true":
        init_sample_data()
    from app.database import SessionLocal
    db = SessionLocal()
    try:
        seed_default_rules(db)
        list_reference_signatures(db)
    finally:
        db.close()


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
