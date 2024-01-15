"""
Per-position diversity statistics over a set of equal-length (aligned)
sequences: Shannon entropy and a simple majority consensus.

Used by the phylogeny app's "mutation hotspot" view and to derive a
reference/consensus sequence for a clade.
"""

from __future__ import annotations

import math
from collections import Counter
from typing import Sequence as Seq


def shannon_entropy(column: Seq[str]) -> float:
    """Shannon entropy (in bits) of a single alignment column.

    Gaps and ambiguous characters are ignored. A perfectly conserved column
    has entropy 0; an even mix of all four nucleotides has entropy 2.
    """
    counts = Counter(c.upper() for c in column if c.upper() in "ACGTU")
    total = sum(counts.values())
    if total == 0:
        return 0.0
    h = 0.0
    for n in counts.values():
        p = n / total
        h -= p * math.log2(p)
    return h


def entropy_profile(sequences: Seq[str]) -> list[float]:
    """Entropy at each position across equal-length sequences."""
    if not sequences:
        return []
    length = len(sequences[0])
    if any(len(s) != length for s in sequences):
        raise ValueError("entropy_profile requires equal-length sequences")
    profile = []
    for pos in range(length):
        column = [s[pos] for s in sequences]
        profile.append(shannon_entropy(column))
    return profile


def consensus(sequences: Seq[str]) -> str:
    """Majority-rule consensus over equal-length sequences."""
    if not sequences:
        return ""
    length = len(sequences[0])
    out = []
    for pos in range(length):
        counts = Counter(s[pos].upper() for s in sequences if s[pos] != "-")
        out.append(counts.most_common(1)[0][0] if counts else "-")
    return "".join(out)


def hotspots(profile: list[float], threshold: float = 1.0) -> list[int]:
    """Indices whose entropy exceeds ``threshold`` bits — variable sites."""
    return [i for i, h in enumerate(profile) if h >= threshold]
