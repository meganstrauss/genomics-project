"""
Minimal SQLite helper shared by both apps.

Uses thread-local connections plus a process-wide write lock so that an HTTP
request thread and a background job thread never issue interleaved writes
(which can corrupt WAL state). Reads may proceed without the lock.
"""

from __future__ import annotations

import sqlite3
import threading
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator


class Database:
    def __init__(self, path: str | Path, schema: str):
        self.path = str(path)
        self._local = threading.local()
        self._write_lock = threading.RLock()
        Path(self.path).parent.mkdir(parents=True, exist_ok=True)
        with self.write() as conn:
            conn.executescript(schema)

    def _conn(self) -> sqlite3.Connection:
        conn = getattr(self._local, "conn", None)
        if conn is None:
            conn = sqlite3.connect(self.path, check_same_thread=False, timeout=30.0)
            conn.row_factory = sqlite3.Row
            conn.execute("PRAGMA foreign_keys = ON")
            conn.execute("PRAGMA journal_mode = WAL")
            conn.execute("PRAGMA busy_timeout = 30000")
            self._local.conn = conn
        return conn

    @contextmanager
    def write(self) -> Iterator[sqlite3.Connection]:
        conn = self._conn()
        self._write_lock.acquire()
        try:
            yield conn
            conn.commit()
        except Exception:
            conn.rollback()
            raise
        finally:
            self._write_lock.release()

    @contextmanager
    def read(self) -> Iterator[sqlite3.Connection]:
        yield self._conn()
