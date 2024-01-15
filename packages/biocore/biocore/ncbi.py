"""
NCBI Entrez E-utilities client.

Implements the real esearch -> efetch workflow against NCBI's public API
(https://www.ncbi.nlm.nih.gov/books/NBK25501/), including the politeness
rate limit (3 requests/sec without an API key, 10/sec with one).

Because automated bulk access requires network egress that is often
unavailable (CI, sandboxes, offline dev), the client degrades to a
:class:`SyntheticEntrez` source that produces reproducible, biologically
plausible records. Both satisfy the same :class:`EntrezSource` protocol, so
nothing downstream knows or cares which is in use.
"""

from __future__ import annotations

import os
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from typing import Optional, Protocol

from .sequences import Sequence, parse_fasta

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@dataclass
class Query:
    term: str                       # Entrez query, e.g. 'Influenza A virus[Organism]'
    db: str = "nuccore"             # nucleotide database
    retmax: int = 50                # max records
    organism: Optional[str] = None


class EntrezSource(Protocol):
    def search_and_fetch(self, query: Query) -> list[Sequence]: ...


class NcbiEntrez:
    """Live client. Honours NCBI rate limits and optional API key."""

    def __init__(self, api_key: Optional[str] = None, allow_network: bool = False):
        self.api_key = api_key or os.environ.get("NCBI_API_KEY")
        self.allow_network = allow_network
        self._min_interval = 0.34 if not self.api_key else 0.11
        self._last_call = 0.0

    def _throttle(self) -> None:
        elapsed = time.monotonic() - self._last_call
        if elapsed < self._min_interval:
            time.sleep(self._min_interval - elapsed)
        self._last_call = time.monotonic()

    def _get(self, endpoint: str, params: dict) -> str:
        if not self.allow_network:
            raise RuntimeError(
                "NcbiEntrez requires allow_network=True. Use SyntheticEntrez "
                "for offline runs, or set HDT_ALLOW_NETWORK=1 with egress to NCBI."
            )
        if self.api_key:
            params["api_key"] = self.api_key
        self._throttle()
        url = f"{EUTILS_BASE}/{endpoint}?{urllib.parse.urlencode(params)}"
        with urllib.request.urlopen(url, timeout=30) as resp:  # noqa: S310
            return resp.read().decode("utf-8", errors="replace")

    def esearch(self, query: Query) -> list[str]:
        """Return a list of UIDs matching the query term."""
        xml = self._get(
            "esearch.fcgi",
            {"db": query.db, "term": query.term, "retmax": str(query.retmax)},
        )
        # minimal XML scrape to avoid an lxml dependency
        ids = []
        for chunk in xml.split("<Id>")[1:]:
            ids.append(chunk.split("</Id>")[0].strip())
        return ids

    def efetch(self, ids: list[str], db: str = "nuccore") -> list[Sequence]:
        if not ids:
            return []
        fasta = self._get(
            "efetch.fcgi",
            {"db": db, "id": ",".join(ids), "rettype": "fasta", "retmode": "text"},
        )
        return list(parse_fasta(fasta))

    def search_and_fetch(self, query: Query) -> list[Sequence]:
        return self.efetch(self.esearch(query), db=query.db)


def make_source(allow_network: bool = False, synthetic_profile: str = "influenza") -> EntrezSource:
    """Factory: live client when network is allowed, else the synthetic one."""
    if allow_network:
        return NcbiEntrez(allow_network=True)
    from .synthetic import SyntheticEntrez

    return SyntheticEntrez(profile=synthetic_profile)
