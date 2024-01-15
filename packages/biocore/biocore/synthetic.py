"""
Offline synthetic data sources.

These reproduce the *structure* real queries return, so the full pipeline
runs without network access:

  * SyntheticEntrez (influenza profile): generates viral sequences along a
    branching lineage tree that accumulates substitutions over simulated
    years, so neighbor-joining recovers a sensible time-structured tree.
  * SyntheticAMR: generates a catalogue of resistance genes and a set of
    bacterial genomes carrying subsets of them, with deliberate
    co-occurrence structure (genes that travel together on the same mobile
    elements) for the network analysis to find.

All generation is seeded and deterministic.
"""

from __future__ import annotations

import hashlib
import random
from dataclasses import dataclass, field
from typing import Optional

from .ncbi import Query
from .sequences import Sequence

_NT = "ACGT"


def _seeded(*parts: object) -> random.Random:
    key = "::".join(str(p) for p in parts)
    digest = hashlib.sha256(key.encode()).hexdigest()
    return random.Random(int(digest[:16], 16))


# --------------------------------------------------------------------------
# Phylogeny: time-structured viral sequences
# --------------------------------------------------------------------------

# A small set of named clades that split off in successive eras, loosely
# echoing how influenza A H3N2 clades have replaced one another over time.
_CLADES = [
    ("Pre-1990", 1985, 1992),
    ("Sydney-like", 1993, 2000),
    ("Fujian-like", 2001, 2008),
    ("Perth-like", 2009, 2014),
    ("HongKong-like", 2015, 2022),
]

_COUNTRIES = ["USA", "China", "Australia", "United Kingdom", "Japan", "Brazil"]
_SEQ_LEN = 900  # ~ a single gene segment (HA-ish), kept short for speed


class SyntheticEntrez:
    def __init__(self, profile: str = "influenza", seed: int = 1933):
        # 1933: first human influenza A virus isolated.
        self.profile = profile
        self.seed = seed
        self._ancestor = "".join(
            _seeded(seed, "ancestor").choice(_NT) for _ in range(_SEQ_LEN)
        )

    def _clade_sequence(self, clade: str) -> str:
        """Each clade descends from the ancestor with a fixed set of clade-
        defining substitutions, so members of a clade cluster together."""
        rng = _seeded(self.seed, "clade", clade)
        bases = list(self._ancestor)
        for _ in range(40):  # clade-defining mutations
            pos = rng.randrange(_SEQ_LEN)
            bases[pos] = rng.choice(_NT)
        return "".join(bases)

    def search_and_fetch(self, query: Query) -> list[Sequence]:
        n = query.retmax
        out: list[Sequence] = []
        for i in range(n):
            rng = _seeded(self.seed, query.term, i)
            clade_name, y0, y1 = rng.choice(_CLADES)
            year = rng.randint(y0, y1)
            clade_seq = list(self._clade_sequence(clade_name))
            # within-clade drift proportional to position in the clade's era
            drift = (year - y0) / max(1, (y1 - y0))
            n_mut = int(rng.uniform(2, 10) * (0.5 + drift))
            for _ in range(n_mut):
                pos = rng.randrange(_SEQ_LEN)
                clade_seq[pos] = rng.choice(_NT)
            country = rng.choice(_COUNTRIES)
            acc = f"SYN{i:05d}"
            out.append(
                Sequence(
                    accession=acc,
                    seq="".join(clade_seq),
                    organism="Influenza A virus (synthetic)",
                    country=country,
                    year=year,
                    metadata={"clade": clade_name},
                )
            )
        return out


# --------------------------------------------------------------------------
# AMR: resistance genes and the genomes that carry them
# --------------------------------------------------------------------------


@dataclass
class AMRGene:
    name: str
    drug_class: str
    seq: str


@dataclass
class BacterialGenome:
    accession: str
    organism: str
    country: str
    year: int
    gene_names: list[str] = field(default_factory=list)


# Real AMR gene families and the drug classes they confer resistance to.
# (Names are real gene families; sequences are synthetic.)
_GENE_CATALOGUE = [
    ("blaTEM", "beta-lactam"),
    ("blaCTX-M", "beta-lactam"),
    ("blaKPC", "carbapenem"),
    ("blaNDM", "carbapenem"),
    ("mecA", "beta-lactam"),
    ("vanA", "glycopeptide"),
    ("vanB", "glycopeptide"),
    ("tetA", "tetracycline"),
    ("tetM", "tetracycline"),
    ("sul1", "sulfonamide"),
    ("sul2", "sulfonamide"),
    ("aac(6')", "aminoglycoside"),
    ("aph(3')", "aminoglycoside"),
    ("qnrS", "fluoroquinolone"),
    ("ermB", "macrolide"),
    ("mcr-1", "polymyxin"),
]

# Mobile-element "cassettes": genes that physically travel together and so
# co-occur far more often than chance. The analysis should recover these.
_CASSETTES = [
    ["blaCTX-M", "sul1", "tetA", "aac(6')"],     # classic ESBL plasmid
    ["blaKPC", "blaNDM", "qnrS"],                 # carbapenemase cluster
    ["mecA", "ermB", "tetM"],                     # MRSA-associated
    ["vanA", "vanB"],                             # vancomycin resistance operon
    ["sul2", "aph(3')", "tetA"],                  # broad plasmid
]

_ORGANISMS = [
    "Escherichia coli",
    "Klebsiella pneumoniae",
    "Staphylococcus aureus",
    "Enterococcus faecium",
    "Pseudomonas aeruginosa",
    "Acinetobacter baumannii",
]


class SyntheticAMR:
    def __init__(self, seed: int = 1928):
        # 1928: Fleming discovers penicillin.
        self.seed = seed

    def genes(self) -> list[AMRGene]:
        out = []
        for name, drug in _GENE_CATALOGUE:
            rng = _seeded(self.seed, "gene", name)
            length = rng.randint(600, 1200)
            seq = "".join(rng.choice(_NT) for _ in range(length))
            out.append(AMRGene(name=name, drug_class=drug, seq=seq))
        return out

    def genomes(self, n: int = 120) -> list[BacterialGenome]:
        out = []
        for i in range(n):
            rng = _seeded(self.seed, "genome", i)
            organism = rng.choice(_ORGANISMS)
            country = rng.choice(_COUNTRIES)
            year = rng.randint(2008, 2023)
            carried: set[str] = set()
            # each genome picks 1-2 cassettes plus a little background noise
            for cassette in rng.sample(_CASSETTES, k=rng.randint(1, 2)):
                # carry most (not always all) of a cassette's genes
                for g in cassette:
                    if rng.random() < 0.85:
                        carried.add(g)
            for name, _ in _GENE_CATALOGUE:
                if rng.random() < 0.05:  # sporadic background carriage
                    carried.add(name)
            out.append(
                BacterialGenome(
                    accession=f"GCA_{900000 + i}",
                    organism=organism,
                    country=country,
                    year=year,
                    gene_names=sorted(carried),
                )
            )
        return out
