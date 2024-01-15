"""Phylogeny Builder HTTP API."""

from __future__ import annotations

import os
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from biocore import Database, JobQueue, JobStatus

from .service import PhylogenyService, SCHEMA

DB_PATH = os.environ.get("PHYLO_DB_PATH", "data/phylo.db")
ALLOW_NETWORK = os.environ.get("PGS_ALLOW_NETWORK", "0") == "1"
K = int(os.environ.get("PHYLO_K", "12"))

app = FastAPI(title="Pathogen Phylogeny Builder", version="0.1.0")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

_db = Database(DB_PATH, SCHEMA)
_queue = JobQueue(max_workers=2)
_svc = PhylogenyService(_db, _queue, allow_network=ALLOW_NETWORK, k=K)


@app.get("/api/health")
def health() -> dict:
    return {"status": "ok", "network": ALLOW_NETWORK}


@app.post("/api/ingest")
def ingest(term: str = "Influenza A virus[Organism]", retmax: int = 40,
           query_tag: str = "default") -> dict:
    r = _svc.ingest(term, retmax, query_tag)
    return r.__dict__


@app.post("/api/build")
def build(query_tag: str = "default") -> dict:
    job_id = _svc.submit_build(query_tag)
    return {"job_id": job_id, "status": "submitted"}


@app.get("/api/jobs/{job_id}")
def job(job_id: str) -> dict:
    j = _queue.get(job_id)
    if not j:
        raise HTTPException(404, "job not found")
    return {"id": j.id, "name": j.name, "status": j.status.value,
            "progress": j.progress, "error": j.error,
            "result": j.result if j.status == JobStatus.SUCCEEDED else None}


@app.get("/api/tree")
def tree(query_tag: str = "default") -> dict:
    t = _svc.latest_tree(query_tag)
    if not t:
        raise HTTPException(404, "no tree built yet — run ingest then build")
    return t


@app.get("/api/timeline")
def timeline(query_tag: str = "default") -> dict:
    return _svc.timeline(query_tag)


_WEB = Path(__file__).resolve().parent / "web"
if _WEB.exists():
    @app.get("/")
    def index() -> FileResponse:
        return FileResponse(_WEB / "index.html")

    app.mount("/", StaticFiles(directory=str(_WEB)), name="web")
