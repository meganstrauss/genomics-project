"""AMR Gene Atlas HTTP API."""

from __future__ import annotations

import os
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from biocore import Database, JobQueue, JobStatus

from .service import AMRService, SCHEMA

DB_PATH = os.environ.get("AMR_DB_PATH", "data/amr.db")
K = int(os.environ.get("AMR_K", "8"))

app = FastAPI(title="AMR Gene Atlas", version="0.1.0")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

_db = Database(DB_PATH, SCHEMA)
_queue = JobQueue(max_workers=2)
_svc = AMRService(_db, _queue, k=K)


@app.get("/api/health")
def health() -> dict:
    return {"status": "ok"}


@app.post("/api/pipeline")
def pipeline(n_genomes: int = 120) -> dict:
    job_id = _svc.submit_pipeline(n_genomes)
    return {"job_id": job_id, "status": "submitted"}


@app.get("/api/jobs/{job_id}")
def job(job_id: str) -> dict:
    j = _queue.get(job_id)
    if not j:
        raise HTTPException(404, "job not found")
    return {"id": j.id, "name": j.name, "status": j.status.value,
            "progress": j.progress, "error": j.error,
            "result": j.result if j.status == JobStatus.SUCCEEDED else None}


@app.get("/api/network")
def network(min_lift: float = 1.5, min_count: int = 5) -> dict:
    return _svc.network(min_lift=min_lift, min_count=min_count)


@app.get("/api/summary")
def summary() -> dict:
    return _svc.summary()


_WEB = Path(__file__).resolve().parent / "web"
if _WEB.exists():
    @app.get("/")
    def index() -> FileResponse:
        return FileResponse(_WEB / "index.html")

    app.mount("/", StaticFiles(directory=str(_WEB)), name="web")
