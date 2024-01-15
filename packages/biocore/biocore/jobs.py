"""In-process background job queue (thread pool + status registry)."""

from __future__ import annotations

import threading
import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Optional


class JobStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"


@dataclass
class Job:
    id: str
    name: str
    status: JobStatus = JobStatus.PENDING
    result: Any = None
    error: Optional[str] = None
    progress: dict = field(default_factory=dict)


class JobQueue:
    def __init__(self, max_workers: int = 2):
        self._executor = ThreadPoolExecutor(max_workers=max_workers)
        self._jobs: dict[str, Job] = {}
        self._lock = threading.Lock()

    def submit(self, name: str, fn: Callable[[Job], Any]) -> str:
        job = Job(id=uuid.uuid4().hex[:12], name=name)
        with self._lock:
            self._jobs[job.id] = job

        def _run() -> None:
            job.status = JobStatus.RUNNING
            try:
                job.result = fn(job)
                job.status = JobStatus.SUCCEEDED
            except Exception as exc:  # noqa: BLE001
                job.status = JobStatus.FAILED
                job.error = f"{exc}\n{traceback.format_exc()}"

        self._executor.submit(_run)
        return job.id

    def get(self, job_id: str) -> Optional[Job]:
        with self._lock:
            return self._jobs.get(job_id)

    def all(self) -> list[Job]:
        with self._lock:
            return list(self._jobs.values())
