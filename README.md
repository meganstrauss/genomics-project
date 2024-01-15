# Pathogen Genomics Suite

Two interactive bioinformatics tools built on a shared analysis core, designed
to pull from public sequence databases (NCBI Entrez / GenBank) and turn raw
genomic data into something you can explore in the browser.

- **Phylogeny Builder** — reconstructs viral lineage trees from genomic
  distance and replays their divergence across collection years.
- **AMR Gene Atlas** — maps which antibiotic-resistance genes co-occur across
  bacterial genomes, surfacing the gene cassettes that drive multidrug
  resistance.

Both run **fully offline** out of the box using a deterministic synthetic data
generator, and both can be pointed at live NCBI data by flipping one flag. No
external map tiles, no paid APIs, no build step for the frontends.

## Why this exists

I wanted to create a visualization tool for pulling from a large datasource. This tool applies analysis directly to the data and lets you visualize it all in one package.

## Architecture

```
pathogen-genomics-suite/
├── packages/biocore/          # shared, installable analysis library
│   └── biocore/
│       ├── sequences.py        # Sequence model + tolerant FASTA parser
│       ├── distance.py         # k-mer Jaccard, p-distance, Jukes-Cantor
│       ├── phylo.py            # neighbor-joining + Newick (from scratch)
│       ├── entropy.py          # Shannon entropy, consensus, hotspots
│       ├── ncbi.py             # NCBI E-utilities client (esearch/efetch)
│       ├── synthetic.py        # offline data generators
│       ├── db.py               # SQLite helper (thread-safe writes)
│       └── jobs.py             # background job queue
├── apps/
│   ├── phylogeny/              # Phylogeny Builder (service + API + web)
│   └── amr_atlas/              # AMR Gene Atlas (service + API + web)
├── tests/                      # pytest suite
├── cli.py                      # headless runner + server launcher
├── Dockerfile / docker-compose.yml
└── .github/workflows/ci.yml
```

## Quick start

```bash
# 1. install the shared library and dependencies
pip install -e packages/biocore
pip install -r requirements.txt

# 2. run a pipeline headless to see it work
python cli.py phylo --retmax 40 --newick
python cli.py amr   --genomes 150 --top 10

# 3. or launch the interactive apps (separate terminals)
python cli.py serve phylo --port 8000   # http://127.0.0.1:8000
python cli.py serve amr   --port 8001   # http://127.0.0.1:8001
```

With Docker:

```bash
docker compose up --build
# phylogeny -> http://localhost:8000
# amr_atlas -> http://localhost:8001
```

In each web app, click **Run pipeline** to ingest data, build the
tree / network, and render it.


## The Phylogeny Builder

**Pipeline:** `esearch`/`efetch` sequences → k-mer Jaccard distance matrix →
neighbor-joining tree → radial time-tree.

- Distances are **alignment-free** (k-mer Jaccard, Mash-style), so no multiple
  sequence alignment is required and indels don't break the comparison.
- The tree is built with a from-scratch implementation of **neighbor-joining**
  (Saitou & Nei, 1987), verified against the canonical textbook example
  (see `tests/`).
- Leaves are keyed by collection year, so the **timeline slider** reveals
  lineages in the order they actually emerged.

The synthetic profile models influenza-like clade replacement: sequences
descend from clade-defining mutations and drift within their era, so NJ
recovers clean, time-ordered clades.

## The AMR Gene Atlas

**Pipeline:** load resistance-gene catalogue + bacterial genomes → single-
linkage gene clustering by sequence similarity → co-occurrence network.

- Each network edge is weighted by **lift** — observed co-occurrence divided by
  what you'd expect if two genes were carried independently. This isolates genes
  that physically travel together (same plasmid / mobile element) from genes
  that are merely individually common.
- The synthetic profile plants real-world cassettes (an ESBL plasmid, a
  carbapenemase cluster, an MRSA-associated set) so the network has genuine
  structure to recover. The carbapenemase cluster (blaKPC / blaNDM / qnrS)
  reliably tops the lift ranking.

## Using live NCBI data

The phylogeny app can pull real sequences from NCBI Entrez:

```bash
export PGS_ALLOW_NETWORK=1
export NCBI_API_KEY=your_key   # optional; raises the rate limit to 10 req/s
python cli.py serve phylo --port 8000
```

The client (`biocore/ncbi.py`) implements the standard `esearch` → `efetch`
workflow and honours NCBI's politeness rate limits. Without the flag it falls
back to the synthetic source, so CI and offline development always work.

## Testing

```bash
pytest tests/ -v
```

The suite covers the parts where correctness matters: distance metrics, the
neighbor-joining branch lengths (against the textbook example), Shannon entropy,
and FASTA parsing.