"""
Sequence distance metrics and distance-matrix construction.

Two families are provided:

  * Alignment-free (k-mer Jaccard / Mash-style): fast, robust to indels,
    used as the default for tree building because it does not require a
    multiple sequence alignment.
  * Alignment-based (p-distance / Jukes-Cantor): the classic proportion of
    differing sites, optionally corrected for multiple substitutions, used
    when sequences are already the same length (e.g. an alignment column
    block).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import combinations
from typing import Callable, Sequence as Seq

from .sequences import Sequence


def kmer_set(seq: str, k: int) -> frozenset[str]:
    s = seq.upper()
    return frozenset(s[i : i + k] for i in range(len(s) - k + 1))


def jaccard_distance(a: frozenset[str], b: frozenset[str]) -> float:
    if not a and not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return 1.0 - (inter / union) if union else 1.0


def p_distance(a: str, b: str) -> float:
    """Proportion of sites that differ, over aligned positions where both
    sequences have a non-gap character. Requires equal length.
    """
    if len(a) != len(b):
        raise ValueError("p-distance requires aligned (equal-length) sequences")
    diffs = 0
    valid = 0
    for x, y in zip(a.upper(), b.upper()):
        if x in "-." or y in "-.":
            continue
        valid += 1
        if x != y:
            diffs += 1
    return diffs / valid if valid else 0.0


def jukes_cantor(p: float) -> float:
    """Jukes-Cantor correction of a p-distance to expected substitutions/site.

    Diverges as p -> 0.75; we cap just below to keep the tree finite.
    """
    p = min(p, 0.7499)
    return -0.75 * math.log(1.0 - (4.0 / 3.0) * p)


@dataclass
class DistanceMatrix:
    labels: list[str]
    matrix: list[list[float]]

    def __post_init__(self) -> None:
        n = len(self.labels)
        if len(self.matrix) != n or any(len(row) != n for row in self.matrix):
            raise ValueError("matrix shape does not match number of labels")

    @property
    def n(self) -> int:
        return len(self.labels)

    def get(self, a: str, b: str) -> float:
        i, j = self.labels.index(a), self.labels.index(b)
        return self.matrix[i][j]

    def to_dict(self) -> dict:
        return {"labels": self.labels, "matrix": self.matrix}


def kmer_distance_matrix(sequences: Seq[Sequence], k: int = 12) -> DistanceMatrix:
    """Build an all-pairs k-mer Jaccard distance matrix.

    k-mer sets are computed once per sequence (O(n) work) and then compared
    pairwise (O(n^2) set intersections), which is far cheaper than O(n^2)
    pairwise alignment for hundreds of sequences.
    """
    labels = [s.accession for s in sequences]
    profiles = [kmer_set(s.cleaned(), k) for s in sequences]
    n = len(labels)
    matrix = [[0.0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        d = jaccard_distance(profiles[i], profiles[j])
        matrix[i][j] = matrix[j][i] = d
    return DistanceMatrix(labels, matrix)


def alignment_distance_matrix(
    sequences: Seq[Sequence], correct: bool = True
) -> DistanceMatrix:
    """Distance matrix from already-aligned, equal-length sequences."""
    labels = [s.accession for s in sequences]
    seqs = [s.seq for s in sequences]
    n = len(labels)
    matrix = [[0.0] * n for _ in range(n)]
    transform: Callable[[float], float] = jukes_cantor if correct else (lambda p: p)
    for i, j in combinations(range(n), 2):
        p = p_distance(seqs[i], seqs[j])
        d = transform(p)
        matrix[i][j] = matrix[j][i] = d
    return DistanceMatrix(labels, matrix)
