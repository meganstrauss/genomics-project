"""
Sequence types and FASTA I/O.

A :class:`Sequence` is the unit of analysis across both apps: a nucleotide or
protein string plus the metadata we can recover from a GenBank/NCBI record
(accession, organism, collection date, country, gene). The FASTA parser is
deliberately tolerant of the messy headers real records carry.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterator, Optional

_NT = frozenset("ACGTU")
_GAP = frozenset("-.")


@dataclass
class Sequence:
    accession: str
    seq: str
    organism: Optional[str] = None
    gene: Optional[str] = None
    country: Optional[str] = None
    year: Optional[int] = None
    metadata: dict = field(default_factory=dict)

    @property
    def length(self) -> int:
        return len(self.seq)

    def cleaned(self) -> str:
        """Uppercase, drop gaps and anything outside the nucleotide alphabet."""
        return "".join(c for c in self.seq.upper() if c in _NT)

    def is_nucleotide(self) -> bool:
        sample = self.seq.upper()[:200]
        letters = [c for c in sample if c.isalpha()]
        if not letters:
            return True
        nt = sum(1 for c in letters if c in _NT)
        return nt / len(letters) > 0.85


def parse_fasta(text: str) -> Iterator[Sequence]:
    """Yield :class:`Sequence` records from FASTA text.

    Header parsing tries the common NCBI conventions but never fails on an
    unrecognised header — it just leaves the metadata fields as ``None``.
    """
    header: Optional[str] = None
    chunks: list[str] = []

    def emit(h: str, body: str) -> Sequence:
        return _sequence_from_header(h, body)

    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                yield emit(header, "".join(chunks))
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line.strip())
    if header is not None:
        yield emit(header, "".join(chunks))


def to_fasta(sequences: list[Sequence], width: int = 70) -> str:
    """Serialise sequences back to FASTA, wrapping body lines at ``width``."""
    out: list[str] = []
    for s in sequences:
        meta_bits = [s.accession]
        if s.organism:
            meta_bits.append(s.organism)
        if s.country:
            meta_bits.append(f"country={s.country}")
        if s.year:
            meta_bits.append(f"year={s.year}")
        out.append(">" + " | ".join(meta_bits))
        body = s.seq
        for i in range(0, len(body), width):
            out.append(body[i : i + width])
    return "\n".join(out) + "\n"


# NCBI headers vary, but year and country are frequently embedded as
# "A/Sydney/5/1997" style strain names or as "/country=..." qualifiers. We
# extract a 4-digit year and a trailing token where we can.
_YEAR_RE = re.compile(r"(?:^|[^\d])((?:19|20)\d{2})(?:[^\d]|$)")
_COUNTRY_RE = re.compile(r"country[=:]\s*([A-Za-z ]+)")


def _sequence_from_header(header: str, body: str) -> Sequence:
    accession = header.split()[0].split("|")[0] if header else "UNKNOWN"
    year = None
    m = _YEAR_RE.search(header)
    if m:
        year = int(m.group(1))
    country = None
    mc = _COUNTRY_RE.search(header)
    if mc:
        country = mc.group(1).strip()
    organism = None
    if "|" in header:
        parts = [p.strip() for p in header.split("|")]
        if len(parts) > 1:
            organism = parts[1]
    return Sequence(
        accession=accession,
        seq=body,
        organism=organism,
        country=country,
        year=year,
    )
