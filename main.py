import os
import uuid
from contextlib import asynccontextmanager

from fastapi import FastAPI, BackgroundTasks, Depends, HTTPException
from fastapi.responses import JSONResponse, Response
from pydantic import BaseModel
from sqlalchemy.ext.asyncio import AsyncSession

from database import init_db, AsyncSessionLocal
from database.models import PipelineJob
from pipeline_controller import run_pipeline
from report_generator.controller import generate_json_report, generate_pdf_report


# ---------------------------------------------------------------------------
# Lifespan – create DB tables on startup
# ---------------------------------------------------------------------------
@asynccontextmanager
async def lifespan(app: FastAPI):
    await init_db()
    yield


app = FastAPI(
    title="Multi-Epitope Vaccine (MEV) Design Pipeline",
    description=(
        "12-step MEV design pipeline: sequence retrieval → antigenicity gate → "
        "physicochemical analysis → epitope prediction → population coverage → "
        "safety filtering → IFN-gamma → docking → model evaluation → "
        "construct design → codon optimization → in-silico cloning."
    ),
    version="2.0.0",
    lifespan=lifespan,
)


# ---------------------------------------------------------------------------
# Dependency – DB session per request
# ---------------------------------------------------------------------------
async def get_db():
    async with AsyncSessionLocal() as session:
        yield session


# ---------------------------------------------------------------------------
# Schemas
# ---------------------------------------------------------------------------
class PipelineRunRequest(BaseModel):
    uniprot_id: str
    host_organism: str = "E. coli"


class PipelineRunResponse(BaseModel):
    job_id: str
    status: str
    message: str


class JobStatusResponse(BaseModel):
    job_id: str
    uniprot_id: str
    status: str
    created_at: str | None = None
    completed_at: str | None = None
    error_message: str | None = None


# ---------------------------------------------------------------------------
# Background task wrapper
# ---------------------------------------------------------------------------
async def _run_pipeline_bg(job_id: str, uniprot_id: str):
    """Run the pipeline in a background context with its own DB session."""
    async with AsyncSessionLocal() as db:
        await run_pipeline(job_id, uniprot_id, db)


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------
@app.get("/api/health")
async def health_check():
    return {"status": "healthy", "service": "immunoinformatics-pipeline"}


@app.post("/api/pipeline/run", response_model=PipelineRunResponse)
async def start_pipeline(
    request: PipelineRunRequest,
    background_tasks: BackgroundTasks,
    db: AsyncSession = Depends(get_db),
):
    """Start a new pipeline evaluation for the given UniProt ID."""
    job_id = str(uuid.uuid4())

    job = PipelineJob(
        id=job_id,
        uniprot_id=request.uniprot_id.strip().upper(),
        status="pending",
    )
    db.add(job)
    await db.commit()

    background_tasks.add_task(_run_pipeline_bg, job_id, request.uniprot_id.strip().upper())

    return PipelineRunResponse(
        job_id=job_id,
        status="pending",
        message="Pipeline started. Poll /api/pipeline/{job_id}/status for updates.",
    )


@app.get("/api/pipeline/{job_id}/status", response_model=JobStatusResponse)
async def get_job_status(job_id: str, db: AsyncSession = Depends(get_db)):
    """Poll the current status of a pipeline job."""
    job = await db.get(PipelineJob, job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found.")

    return JobStatusResponse(
        job_id=job.id,
        uniprot_id=job.uniprot_id,
        status=job.status,
        created_at=job.created_at.isoformat() if job.created_at else None,
        completed_at=job.completed_at.isoformat() if job.completed_at else None,
        error_message=job.error_message,
    )


@app.get("/api/pipeline/{job_id}/report")
async def get_json_report(job_id: str, db: AsyncSession = Depends(get_db)):
    """Return the full structured JSON report for a completed job."""
    job = await db.get(PipelineJob, job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found.")
    if job.status != "completed":
        raise HTTPException(
            status_code=400,
            detail=f"Job is not yet completed. Current status: {job.status}",
        )

    report = await generate_json_report(job_id, db)
    return JSONResponse(content=report)


@app.get("/api/pipeline/{job_id}/report/pdf")
async def get_pdf_report(job_id: str, db: AsyncSession = Depends(get_db)):
    """Download the PDF vaccine candidacy report."""
    job = await db.get(PipelineJob, job_id)
    if not job:
        raise HTTPException(status_code=404, detail="Job not found.")
    if job.status != "completed":
        raise HTTPException(
            status_code=400,
            detail=f"Job is not yet completed. Current status: {job.status}",
        )

    pdf_bytes = await generate_pdf_report(job_id, db)
    return Response(
        content=bytes(pdf_bytes),
        media_type="application/pdf",
        headers={
            "Content-Disposition": f'attachment; filename="report_{job_id}.pdf"'
        },
    )
