"""
Phylogeny Builder service.

Workflow: ingest sequences for a query (organism/year window) -> persist ->
compute a k-mer distance matrix -> build a neighbor-joining tree -> expose
the tree plus a year-keyed timeline so the frontend can animate lineage
divergence over time.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Optional

from biocore import (
    DistanceMatrix,
    Database,
    Job,
    JobQueue,
    Query,
    Sequence,
    kmer_distance_matrix,
    neighbor_joining,
    to_dict,
    to_newick,
)
from biocore.ncbi import EntrezSource, make_source

SCHEMA = """
CREATE TABLE IF NOT EXISTS sequences (
    accession   TEXT PRIMARY KEY,
    organism    TEXT,
    country     TEXT,
    year        INTEGER,
    clade       TEXT,
    length      INTEGER NOT NULL,
    seq         TEXT NOT NULL,
    query_tag   TEXT
);
CREATE INDEX IF NOT EXISTS idx_seq_year ON sequences(year);
CREATE INDEX IF NOT EXISTS idx_seq_tag ON sequences(query_tag);

CREATE TABLE IF NOT EXISTS trees (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    query_tag   TEXT NOT NULL,
    k           INTEGER NOT NULL,
    n_taxa      INTEGER NOT NULL,
    newick      TEXT NOT NULL,
    tree_json   TEXT NOT NULL,
    created_at  TEXT NOT NULL DEFAULT (datetime('now'))
);
"""


@dataclass
class IngestResult:
    query_tag: str
    added: int
    total: int


class PhylogenyService:
    def __init__(self, db: Database, queue: JobQueue, source: Optional[EntrezSource] = None,
                 allow_network: bool = False, k: int = 12):
        self.db = db
        self.queue = queue
        self.k = k
        self.source = source or make_source(allow_network, synthetic_profile="influenza")

    # ---- ingestion ----
    def ingest(self, term: str, retmax: int = 40, query_tag: str = "default") -> IngestResult:
        seqs = self.source.search_and_fetch(Query(term=term, retmax=retmax))
        added = 0
        with self.db.write() as conn:
            for s in seqs:
                cur = conn.execute(
                    """INSERT OR IGNORE INTO sequences
                       (accession, organism, country, year, clade, length, seq, query_tag)
                       VALUES (?,?,?,?,?,?,?,?)""",
                    (s.accession, s.organism, s.country, s.year,
                     s.metadata.get("clade"), s.length, s.seq, query_tag),
                )
                added += cur.rowcount
        with self.db.read() as conn:
            total = conn.execute(
                "SELECT COUNT(*) AS n FROM sequences WHERE query_tag=?", (query_tag,)
            ).fetchone()["n"]
        return IngestResult(query_tag=query_tag, added=added, total=total)

    def _load_sequences(self, query_tag: str) -> list[Sequence]:
        with self.db.read() as conn:
            rows = conn.execute(
                "SELECT * FROM sequences WHERE query_tag=? ORDER BY year, accession",
                (query_tag,),
            ).fetchall()
        return [
            Sequence(
                accession=r["accession"], seq=r["seq"], organism=r["organism"],
                country=r["country"], year=r["year"],
                metadata={"clade": r["clade"]} if r["clade"] else {},
            )
            for r in rows
        ]

    # ---- tree building ----
    def build_tree(self, query_tag: str = "default", job: Optional[Job] = None) -> dict:
        seqs = self._load_sequences(query_tag)
        if len(seqs) < 3:
            raise ValueError("need at least 3 sequences to build a tree")
        if job is not None:
            job.progress = {"stage": "distance matrix", "n": len(seqs)}
        dm: DistanceMatrix = kmer_distance_matrix(seqs, k=self.k)
        if job is not None:
            job.progress = {"stage": "neighbor joining", "n": len(seqs)}
        root = neighbor_joining(dm)
        newick = to_newick(root)
        tree_json = to_dict(root)

        # attach metadata (year, country, clade) to leaves for the renderer
        meta = {s.accession: {"year": s.year, "country": s.country,
                              "clade": s.metadata.get("clade")} for s in seqs}
        _decorate_leaves(tree_json, meta)

        with self.db.write() as conn:
            conn.execute(
                """INSERT INTO trees (query_tag, k, n_taxa, newick, tree_json)
                   VALUES (?,?,?,?,?)""",
                (query_tag, self.k, len(seqs), newick, json.dumps(tree_json)),
            )
        return {"newick": newick, "tree": tree_json, "n_taxa": len(seqs)}

    def submit_build(self, query_tag: str = "default") -> str:
        def _work(job: Job) -> dict:
            return self.build_tree(query_tag, job)
        return self.queue.submit(f"build_tree:{query_tag}", _work)

    def latest_tree(self, query_tag: str = "default") -> Optional[dict]:
        with self.db.read() as conn:
            row = conn.execute(
                "SELECT * FROM trees WHERE query_tag=? ORDER BY id DESC LIMIT 1",
                (query_tag,),
            ).fetchone()
        if not row:
            return None
        return {
            "newick": row["newick"],
            "tree": json.loads(row["tree_json"]),
            "n_taxa": row["n_taxa"],
            "created_at": row["created_at"],
        }

    def timeline(self, query_tag: str = "default") -> dict:
        """Per-year counts and per-clade emergence, for the time animation."""
        with self.db.read() as conn:
            rows = conn.execute(
                "SELECT year, clade, country FROM sequences WHERE query_tag=?",
                (query_tag,),
            ).fetchall()
        by_year: dict[int, int] = {}
        clade_first_year: dict[str, int] = {}
        for r in rows:
            y = r["year"]
            if y is None:
                continue
            by_year[y] = by_year.get(y, 0) + 1
            c = r["clade"]
            if c:
                clade_first_year[c] = min(clade_first_year.get(c, y), y)
        return {
            "by_year": dict(sorted(by_year.items())),
            "clade_emergence": dict(sorted(clade_first_year.items(), key=lambda kv: kv[1])),
        }


def _decorate_leaves(node: dict, meta: dict) -> None:
    if "children" in node:
        for c in node["children"]:
            _decorate_leaves(c, meta)
    elif node.get("name") in meta:
        node.update(meta[node["name"]])
