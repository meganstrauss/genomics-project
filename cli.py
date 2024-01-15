#!/usr/bin/env python3
"""
Pathogen Genomics Suite — command-line interface.

Run either pipeline headless (no web server) for scripting, CI, or quick
inspection. Examples:

    python cli.py phylo --retmax 40 --newick
    python cli.py amr --genomes 150 --top 10
    python cli.py serve phylo --port 8000
    python cli.py serve amr   --port 8001

The pipelines run against synthetic data by default. Set PGS_ALLOW_NETWORK=1
(and have egress to NCBI) to use the live Entrez client for the phylogeny app.
"""

from __future__ import annotations

import argparse
import os
import sys

# make packages importable when run from the repo root
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "packages", "biocore"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "apps"))

from biocore import Database, JobQueue  # noqa: E402


def cmd_phylo(args: argparse.Namespace) -> None:
    from phylogeny.service import PhylogenyService, SCHEMA

    db = Database(args.db, SCHEMA)
    svc = PhylogenyService(db, JobQueue(),
                           allow_network=os.environ.get("PGS_ALLOW_NETWORK") == "1",
                           k=args.k)
    ing = svc.ingest(args.term, retmax=args.retmax, query_tag=args.tag)
    print(f"ingested {ing.added} new sequences ({ing.total} total)")
    tree = svc.build_tree(args.tag)
    print(f"built NJ tree over {tree['n_taxa']} taxa")
    tl = svc.timeline(args.tag)
    print("clade emergence:")
    for clade, year in tl["clade_emergence"].items():
        print(f"  {year}  {clade}")
    if args.newick:
        print("\nNewick:")
        print(tree["newick"])


def cmd_amr(args: argparse.Namespace) -> None:
    from amr_atlas.service import AMRService, SCHEMA

    db = Database(args.db, SCHEMA)
    svc = AMRService(db, JobQueue(), k=args.k)
    ing = svc.ingest(args.genomes)
    print(f"ingested {ing.genes} genes, {ing.genomes} genomes, {ing.presences} presences")
    cl = svc.cluster_genes()
    print(f"clustered into {cl['clusters']} gene clusters")
    net = svc.network()
    print(f"\ntop {args.top} co-occurring cassettes (by lift):")
    for e in net["edges"][: args.top]:
        print(f"  {e['source']:12s} + {e['target']:12s}  lift={e['lift']:>5}  ({e['count']} genomes)")


def cmd_serve(args: argparse.Namespace) -> None:
    import uvicorn

    target = {
        "phylo": "phylogeny.api:app",
        "amr": "amr_atlas.api:app",
    }[args.app]
    uvicorn.run(target, host=args.host, port=args.port, reload=False)


def main() -> None:
    p = argparse.ArgumentParser(description="Pathogen Genomics Suite CLI")
    sub = p.add_subparsers(dest="command", required=True)

    pp = sub.add_parser("phylo", help="run the phylogeny pipeline headless")
    pp.add_argument("--term", default="Influenza A virus[Organism]")
    pp.add_argument("--retmax", type=int, default=40)
    pp.add_argument("--tag", default="cli")
    pp.add_argument("--k", type=int, default=12)
    pp.add_argument("--db", default="data/phylo.db")
    pp.add_argument("--newick", action="store_true", help="print the Newick string")
    pp.set_defaults(func=cmd_phylo)

    pa = sub.add_parser("amr", help="run the AMR pipeline headless")
    pa.add_argument("--genomes", type=int, default=120)
    pa.add_argument("--k", type=int, default=8)
    pa.add_argument("--top", type=int, default=8)
    pa.add_argument("--db", default="data/amr.db")
    pa.set_defaults(func=cmd_amr)

    ps = sub.add_parser("serve", help="serve one of the web apps")
    ps.add_argument("app", choices=["phylo", "amr"])
    ps.add_argument("--host", default="127.0.0.1")
    ps.add_argument("--port", type=int, default=8000)
    ps.set_defaults(func=cmd_serve)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
