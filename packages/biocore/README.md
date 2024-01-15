# biocore

Shared bioinformatics primitives for the Pathogen Genomics Suite. Installable
on its own:

```bash
pip install -e packages/biocore
```

## Modules

| Module        | What it provides                                                        |
| ------------- | ----------------------------------------------------------------------- |
| `sequences`   | `Sequence` dataclass, tolerant FASTA parse/serialize                    |
| `distance`    | k-mer Jaccard, p-distance, Jukes-Cantor, `DistanceMatrix`               |
| `phylo`       | neighbor-joining tree, Newick + JSON serialization                      |
| `entropy`     | per-column Shannon entropy, consensus, variable-site detection          |
| `ncbi`        | NCBI E-utilities client (esearch/efetch) + offline source factory       |
| `synthetic`   | deterministic synthetic data for phylogeny and AMR                      |
| `db`          | thread-safe SQLite helper                                               |
| `jobs`        | in-process background job queue                                         |

## Example

```python
from biocore import Sequence, kmer_distance_matrix, neighbor_joining, to_newick

seqs = [Sequence("a", "ACGT..."), Sequence("b", "ACGA...")]
dm = kmer_distance_matrix(seqs, k=12)
tree = neighbor_joining(dm)
print(to_newick(tree))
```

The neighbor-joining implementation is verified against the standard textbook
example in the suite's test suite.
