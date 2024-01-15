"""biocore — shared bioinformatics primitives for the Pathogen Genomics Suite."""

from .sequences import Sequence, parse_fasta, to_fasta
from .distance import (
    DistanceMatrix,
    kmer_distance_matrix,
    alignment_distance_matrix,
    kmer_set,
    jaccard_distance,
    p_distance,
    jukes_cantor,
)
from .phylo import Node, neighbor_joining, to_newick, to_dict
from .entropy import shannon_entropy, entropy_profile, consensus, hotspots
from .ncbi import Query, NcbiEntrez, make_source

__all__ = [
    "Sequence",
    "parse_fasta",
    "to_fasta",
    "DistanceMatrix",
    "kmer_distance_matrix",
    "alignment_distance_matrix",
    "kmer_set",
    "jaccard_distance",
    "p_distance",
    "jukes_cantor",
    "Node",
    "neighbor_joining",
    "to_newick",
    "to_dict",
    "shannon_entropy",
    "entropy_profile",
    "consensus",
    "hotspots",
    "Query",
    "NcbiEntrez",
    "make_source",
]

from .db import Database  # noqa: E402
__all__.append("Database")

from .jobs import Job, JobQueue, JobStatus  # noqa: E402
__all__ += ["Job", "JobQueue", "JobStatus"]
