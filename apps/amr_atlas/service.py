"""
AMR Gene Atlas service.

Workflow: load the resistance-gene catalogue and a set of bacterial genomes
that carry subsets of those genes -> cluster genes by sequence similarity
(single-linkage over k-mer Jaccard) -> build a co-occurrence network where an
edge between two genes is weighted by how often they appear in the same
genome, relative to chance.

The headline statistic for each edge is a lift / observed-over-expected ratio,
which surfaces genes that physically travel together (e.g. on the same
plasmid) rather than ones that are merely individually common.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from itertools import combinations
from typing import Optional

from biocore import Database, Job, JobQueue, kmer_set, jaccard_distance
from biocore.synthetic import SyntheticAMR

SCHEMA = """
CREATE TABLE IF NOT EXISTS genes (
    name        TEXT PRIMARY KEY,
    drug_class  TEXT NOT NULL,
    length      INTEGER NOT NULL,
    seq         TEXT NOT NULL,
    cluster     INTEGER
);
CREATE TABLE IF NOT EXISTS genomes (
    accession   TEXT PRIMARY KEY,
    organism    TEXT NOT NULL,
    country     TEXT,
    year        INTEGER
);
CREATE TABLE IF NOT EXISTS gene_presence (
    genome      TEXT NOT NULL REFERENCES genomes(accession),
    gene        TEXT NOT NULL REFERENCES genes(name),
    PRIMARY KEY (genome, gene)
);
CREATE INDEX IF NOT EXISTS idx_presence_gene ON gene_presence(gene);
"""


@dataclass
class IngestResult:
    genes: int
    genomes: int
    presences: int


class AMRService:
    def __init__(self, db: Database, queue: JobQueue, k: int = 8,
                 source: Optional[SyntheticAMR] = None):
        self.db = db
        self.queue = queue
        self.k = k
        self.source = source or SyntheticAMR()

    # ---- ingestion ----
    def ingest(self, n_genomes: int = 120) -> IngestResult:
        genes = self.source.genes()
        genomes = self.source.genomes(n_genomes)
        presences = 0
        with self.db.write() as conn:
            for g in genes:
                conn.execute(
                    """INSERT OR IGNORE INTO genes (name, drug_class, length, seq)
                       VALUES (?,?,?,?)""",
                    (g.name, g.drug_class, len(g.seq), g.seq),
                )
            for gn in genomes:
                conn.execute(
                    """INSERT OR IGNORE INTO genomes (accession, organism, country, year)
                       VALUES (?,?,?,?)""",
                    (gn.accession, gn.organism, gn.country, gn.year),
                )
                for gene_name in gn.gene_names:
                    cur = conn.execute(
                        "INSERT OR IGNORE INTO gene_presence (genome, gene) VALUES (?,?)",
                        (gn.accession, gene_name),
                    )
                    presences += cur.rowcount
        return IngestResult(genes=len(genes), genomes=len(genomes), presences=presences)

    # ---- gene clustering (single linkage over sequence similarity) ----
    def cluster_genes(self, threshold: float = 0.6, job: Optional[Job] = None) -> dict:
        with self.db.read() as conn:
            rows = conn.execute("SELECT name, seq FROM genes").fetchall()
        names = [r["name"] for r in rows]
        profiles = {r["name"]: kmer_set(r["seq"].upper(), self.k) for r in rows}

        # union-find for single-linkage clustering
        parent = {n: n for n in names}

        def find(x: str) -> str:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a: str, b: str) -> None:
            parent[find(a)] = find(b)

        for a, b in combinations(names, 2):
            if jaccard_distance(profiles[a], profiles[b]) <= threshold:
                union(a, b)

        # assign integer cluster ids
        roots = {}
        cluster_of = {}
        for n in names:
            r = find(n)
            if r not in roots:
                roots[r] = len(roots)
            cluster_of[n] = roots[r]

        with self.db.write() as conn:
            for n, c in cluster_of.items():
                conn.execute("UPDATE genes SET cluster=? WHERE name=?", (c, n))
        return {"clusters": len(roots), "assignment": cluster_of}

    # ---- co-occurrence network ----
    def network(self, min_lift: float = 1.5, min_count: int = 5) -> dict:
        with self.db.read() as conn:
            gene_rows = conn.execute(
                "SELECT name, drug_class, cluster FROM genes"
            ).fetchall()
            n_genomes = conn.execute("SELECT COUNT(*) AS n FROM genomes").fetchone()["n"]
            presence = conn.execute(
                "SELECT genome, gene FROM gene_presence"
            ).fetchall()

        # per-genome gene sets
        by_genome: dict[str, set[str]] = {}
        gene_count: dict[str, int] = {}
        for r in presence:
            by_genome.setdefault(r["genome"], set()).add(r["gene"])
            gene_count[r["gene"]] = gene_count.get(r["gene"], 0) + 1

        # pairwise co-occurrence
        pair_count: dict[tuple[str, str], int] = {}
        for genes in by_genome.values():
            for a, b in combinations(sorted(genes), 2):
                pair_count[(a, b)] = pair_count.get((a, b), 0) + 1

        nodes = [
            {
                "id": r["name"],
                "drug_class": r["drug_class"],
                "cluster": r["cluster"],
                "prevalence": gene_count.get(r["name"], 0),
            }
            for r in gene_rows
        ]

        edges = []
        N = max(n_genomes, 1)
        for (a, b), c in pair_count.items():
            # expected co-occurrence if independent: P(a)*P(b)*N
            pa = gene_count.get(a, 0) / N
            pb = gene_count.get(b, 0) / N
            expected = pa * pb * N
            lift = (c / expected) if expected > 0 else 0.0
            if c >= min_count and lift >= min_lift:
                edges.append({
                    "source": a, "target": b, "count": c,
                    "lift": round(lift, 3),
                })
        edges.sort(key=lambda e: e["lift"], reverse=True)
        return {"nodes": nodes, "edges": edges, "n_genomes": n_genomes}

    def submit_pipeline(self, n_genomes: int = 120) -> str:
        def _work(job: Job) -> dict:
            job.progress = {"stage": "ingest"}
            ing = self.ingest(n_genomes)
            job.progress = {"stage": "cluster"}
            cl = self.cluster_genes(job=job)
            job.progress = {"stage": "network"}
            net = self.network()
            return {"ingest": ing.__dict__, "clusters": cl["clusters"],
                    "edges": len(net["edges"])}
        return self.queue.submit("amr_pipeline", _work)

    def summary(self) -> dict:
        with self.db.read() as conn:
            genes = conn.execute("SELECT COUNT(*) AS n FROM genes").fetchone()["n"]
            genomes = conn.execute("SELECT COUNT(*) AS n FROM genomes").fetchone()["n"]
            by_class = conn.execute(
                """SELECT g.drug_class, COUNT(*) AS carriage
                   FROM gene_presence p JOIN genes g ON g.name=p.gene
                   GROUP BY g.drug_class ORDER BY carriage DESC"""
            ).fetchall()
        return {
            "genes": genes,
            "genomes": genomes,
            "carriage_by_class": [dict(r) for r in by_class],
        }
