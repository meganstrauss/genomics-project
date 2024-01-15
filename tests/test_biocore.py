"""Tests for the biocore library — the parts where correctness matters."""

import math

import pytest

from biocore import (
    DistanceMatrix,
    Sequence,
    entropy_profile,
    jaccard_distance,
    kmer_set,
    neighbor_joining,
    p_distance,
    parse_fasta,
    shannon_entropy,
    to_newick,
)
from biocore.distance import jukes_cantor, kmer_distance_matrix


# ---- distance metrics ----

def test_jaccard_identical_is_zero():
    a = kmer_set("ACGTACGT", 3)
    assert jaccard_distance(a, a) == 0.0


def test_jaccard_disjoint_is_one():
    a = kmer_set("AAAA", 2)
    b = kmer_set("CCCC", 2)
    assert jaccard_distance(a, b) == 1.0


def test_p_distance_counts_mismatches():
    assert p_distance("ACGT", "ACGA") == pytest.approx(0.25)


def test_p_distance_ignores_gaps():
    assert p_distance("AC-T", "ACGT") == pytest.approx(0.0)


def test_p_distance_requires_equal_length():
    with pytest.raises(ValueError):
        p_distance("ACGT", "ACG")


def test_jukes_cantor_monotonic_and_zero_at_zero():
    assert jukes_cantor(0.0) == pytest.approx(0.0)
    assert jukes_cantor(0.1) < jukes_cantor(0.2)


# ---- neighbor joining: the canonical textbook example ----

def test_neighbor_joining_textbook_branch_lengths():
    labels = ["a", "b", "c", "d", "e"]
    M = [
        [0, 5, 9, 9, 8],
        [5, 0, 10, 10, 9],
        [9, 10, 0, 8, 7],
        [9, 10, 8, 0, 3],
        [8, 9, 7, 3, 0],
    ]
    dm = DistanceMatrix(labels, [[float(x) for x in row] for row in M])
    root = neighbor_joining(dm)
    by_name = {lf.name: lf for lf in root.leaves()}
    # Known result: the first join is (a,b) with branch a=2, b=3.
    assert by_name["a"].branch_length == pytest.approx(2.0, abs=1e-6)
    assert by_name["b"].branch_length == pytest.approx(3.0, abs=1e-6)


def test_neighbor_joining_recovers_all_taxa():
    seqs = [Sequence(f"s{i}", "ACGT" * 10) for i in range(5)]
    dm = kmer_distance_matrix(seqs, k=3)
    root = neighbor_joining(dm)
    assert len(root.leaves()) == 5


def test_newick_is_parenthesised_and_terminated():
    labels = ["x", "y", "z"]
    dm = DistanceMatrix(labels, [[0, 1, 2], [1, 0, 1], [2, 1, 0]])
    nw = to_newick(neighbor_joining(dm))
    assert nw.startswith("(") and nw.endswith(";")
    for label in labels:
        assert label in nw


# ---- entropy ----

def test_shannon_entropy_uniform_is_two_bits():
    assert shannon_entropy(["A", "C", "G", "T"]) == pytest.approx(2.0)


def test_shannon_entropy_conserved_is_zero():
    assert shannon_entropy(["A", "A", "A", "A"]) == pytest.approx(0.0)


def test_entropy_profile_length_matches_alignment():
    profile = entropy_profile(["ACGT", "ACGA", "ACGC"])
    assert len(profile) == 4
    # first three columns conserved, last variable
    assert profile[0] == pytest.approx(0.0)
    assert profile[3] > 0.0


# ---- FASTA ----

def test_parse_fasta_basic():
    text = ">seq1 | Influenza A | year=1999\nACGT\nACGT\n>seq2\nTTTT\n"
    seqs = list(parse_fasta(text))
    assert len(seqs) == 2
    assert seqs[0].accession == "seq1"
    assert seqs[0].seq == "ACGTACGT"
    assert seqs[0].year == 1999


def test_parse_fasta_empty():
    assert list(parse_fasta("")) == []
