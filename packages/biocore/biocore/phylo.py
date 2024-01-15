"""
Neighbor-joining (Saitou & Nei, 1987) phylogenetic tree reconstruction.

Given a pairwise distance matrix, NJ produces an unrooted tree with branch
lengths. The algorithm repeatedly:
  1. computes the Q-matrix Q(i,j) = (n-2)·d(i,j) - Σd(i,·) - Σd(j,·),
  2. joins the pair (i,j) minimising Q into a new node u,
  3. assigns branch lengths from i and j to u,
  4. updates distances from u to every remaining node,
until three nodes remain, which are joined to a final internal node.

We build an explicit tree of :class:`Node` objects and can emit Newick.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from .distance import DistanceMatrix


@dataclass
class Node:
    name: Optional[str] = None          # leaf label, or None for internal
    branch_length: float = 0.0          # length of edge to this node's parent
    children: list["Node"] = field(default_factory=list)
    id: int = 0

    @property
    def is_leaf(self) -> bool:
        return not self.children

    def leaves(self) -> list["Node"]:
        if self.is_leaf:
            return [self]
        out: list[Node] = []
        for c in self.children:
            out.extend(c.leaves())
        return out


def neighbor_joining(dm: DistanceMatrix) -> Node:
    """Construct an NJ tree and return its (arbitrary) root node."""
    if dm.n < 2:
        raise ValueError("need at least two taxa")

    # Working state: active node objects and a mutable distance dict.
    next_id = 0

    def new_id() -> int:
        nonlocal next_id
        next_id += 1
        return next_id

    nodes: list[Node] = [Node(name=label, id=new_id()) for label in dm.labels]
    # distances keyed by node id
    dist: dict[int, dict[int, float]] = {}
    for a in range(dm.n):
        ia = nodes[a].id
        dist[ia] = {}
        for b in range(dm.n):
            dist[ia][nodes[b].id] = dm.matrix[a][b]

    active = list(nodes)

    while len(active) > 2:
        n = len(active)
        # row sums r(i) = Σ_j d(i,j)
        r = {node.id: sum(dist[node.id][o.id] for o in active) for node in active}

        # find pair minimising Q
        best_pair: tuple[Node, Node] | None = None
        best_q = float("inf")
        for x in range(n):
            for y in range(x + 1, n):
                i, j = active[x], active[y]
                q = (n - 2) * dist[i.id][j.id] - r[i.id] - r[j.id]
                if q < best_q:
                    best_q = q
                    best_pair = (i, j)

        assert best_pair is not None
        i, j = best_pair
        dij = dist[i.id][j.id]

        # branch lengths from i and j to the new node u
        delta = (r[i.id] - r[j.id]) / (n - 2)
        li = 0.5 * dij + 0.5 * delta
        lj = dij - li
        # guard against small negative branch lengths from noisy distances
        li = max(li, 0.0)
        lj = max(lj, 0.0)

        u = Node(id=new_id(), children=[i, j])
        i.branch_length = li
        j.branch_length = lj

        # distances from u to each remaining node k
        dist[u.id] = {}
        for k in active:
            if k is i or k is j:
                continue
            duk = 0.5 * (dist[i.id][k.id] + dist[j.id][k.id] - dij)
            dist[u.id][k.id] = duk
            dist[k.id][u.id] = duk
        dist[u.id][u.id] = 0.0

        active = [k for k in active if k is not i and k is not j]
        active.append(u)

    # join the final two nodes
    a, b = active
    dab = dist[a.id][b.id]
    root = Node(id=new_id(), children=[a, b])
    a.branch_length = dab / 2.0
    b.branch_length = dab / 2.0
    return root


def to_newick(root: Node) -> str:
    """Serialise a tree to a Newick string."""

    def fmt(node: Node) -> str:
        if node.is_leaf:
            return f"{node.name}:{node.branch_length:.5f}"
        inner = ",".join(fmt(c) for c in node.children)
        return f"({inner}):{node.branch_length:.5f}"

    # strip the root's own (meaningless) branch length
    body = ",".join(fmt(c) for c in root.children)
    return f"({body});"


def to_dict(root: Node) -> dict:
    """JSON-serialisable nested tree for the frontend renderer."""

    def conv(node: Node) -> dict:
        d = {
            "name": node.name,
            "branch_length": round(node.branch_length, 6),
            "id": node.id,
        }
        if node.children:
            d["children"] = [conv(c) for c in node.children]
        return d

    return conv(root)
