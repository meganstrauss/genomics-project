# Phylogeny Builder

Reconstructs viral lineage trees from genomic distance and animates their
divergence over collection years.

**Pipeline:** NCBI `esearch`/`efetch` → k-mer Jaccard distance matrix →
neighbor-joining tree → radial time-tree with a year slider.

```bash
python cli.py serve phylo --port 8000   # from repo root
```

### API
| Method | Route                    | Purpose                          |
| ------ | ------------------------ | -------------------------------- |
| POST   | `/api/ingest`            | fetch + store sequences          |
| POST   | `/api/build`             | submit NJ tree build (async job) |
| GET    | `/api/jobs/{id}`         | poll a job                       |
| GET    | `/api/tree`              | latest tree (Newick + JSON)      |
| GET    | `/api/timeline`          | per-year counts, clade emergence |
