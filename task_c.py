"""
Part 1 - Task C: Epidemic Models (Troll Propagation)
Models how controversial behaviour spreads through the Wikidata editor network.
"""

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random, os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs_part1")
os.makedirs(OUTPUT_DIR, exist_ok=True)
random.seed(42)
np.random.seed(42)


def propagation_plausibility(G: nx.Graph, node_a: int, node_b: int) -> dict:
    """
    Given two editors (node_a, node_b) who may be 'possibly trolling',
    estimate how plausible it is that the behaviour has NOT yet spread.

    Approach:
      - Shortest path distance: the farther apart, the less likely propagation
      - Number of common neighbours: more shared neighbours → more risk
      - Jaccard similarity: normalised overlap of neighbourhoods
      - Edge weight (shared threads): heavier edges → stronger connection
    """
    Gc = max(nx.connected_components(G), key=len)
    Gc = G.subgraph(Gc).copy()

    if node_a not in Gc or node_b not in Gc:
        nodes = list(Gc.nodes())
        node_a, node_b = nodes[0], nodes[1]

    result = {}

    try:
        path   = nx.shortest_path(Gc, node_a, node_b)
        dist   = len(path) - 1
    except nx.NetworkXNoPath:
        path, dist = [], float("inf")
    result["shortest_path"]     = path
    result["distance"]          = dist

    beta = 0.5
    p_not_propagated = (1 - beta) ** dist if dist < float("inf") else 1.0
    result["p_not_propagated"]  = p_not_propagated

    nbrs_a = set(Gc.neighbors(node_a))
    nbrs_b = set(Gc.neighbors(node_b))
    common = nbrs_a & nbrs_b
    result["common_neighbours"] = len(common)
    result["common_nodes"]      = list(common)[:10]

    union = nbrs_a | nbrs_b
    jaccard = len(common) / len(union) if union else 0.0
    result["jaccard_similarity"] = jaccard

    if Gc.has_edge(node_a, node_b):
        result["direct_edge_weight"] = Gc[node_a][node_b].get("weight", 1)
    else:
        result["direct_edge_weight"] = 0

    print(f"\n  Editor {node_a} ↔ Editor {node_b}")
    print(f"  Shortest path distance : {dist}")
    print(f"  Path                   : {path}")
    print(f"  Common neighbours      : {len(common)}")
    print(f"  Jaccard similarity     : {jaccard:.4f}")
    print(f"  Direct edge weight     : {result['direct_edge_weight']}")
    print(f"  P(not yet propagated)  : {p_not_propagated:.4f}  [β={beta}]")

    return result


def priority_list(G: nx.Graph, infected: list[int], top_n: int = 20) -> pd.DataFrame:
    """
    Given a set of 'infected' (possibly trolling) editors,
    produce a priority list of neighbours to check next.

    Scoring per candidate node u:
      - exposure_score  : sum of edge weights from u to any infected node
      - n_infected_nbrs : count of infected neighbours (cascading risk)
      - degree          : high-degree nodes spread further if infected
      - betweenness     : high-betweenness nodes are structural bridges
      - composite_score : weighted combination
    """
    Gc = max(nx.connected_components(G), key=len)
    Gc = G.subgraph(Gc).copy()

    lcc_nodes = set(Gc.nodes())
    infected  = [n for n in infected if n in lcc_nodes]
    infected_set = set(infected)

    n = Gc.number_of_nodes()
    if n > 3000:
        bc = nx.betweenness_centrality(Gc, normalized=True, k=300)
    else:
        bc = nx.betweenness_centrality(Gc, normalized=True)

    max_bc = max(bc.values()) or 1.0

    candidates = {}
    for inf_node in infected_set:
        for nbr in Gc.neighbors(inf_node):
            if nbr in infected_set:
                continue
            w = Gc[inf_node][nbr].get("weight", 1)
            if nbr not in candidates:
                candidates[nbr] = {"exposure_score": 0, "n_infected_nbrs": 0}
            candidates[nbr]["exposure_score"]  += w
            candidates[nbr]["n_infected_nbrs"] += 1

    rows = []
    for node, scores in candidates.items():
        deg  = Gc.degree(node)
        b    = bc.get(node, 0)
        composite = (
            0.4 * scores["exposure_score"] / max(1, scores["exposure_score"]) +
            0.3 * scores["n_infected_nbrs"] / max(len(infected), 1) +
            0.2 * deg / max(dict(Gc.degree()).values()) +
            0.1 * b / max_bc
        )
        rows.append({
            "node":             node,
            "exposure_score":   scores["exposure_score"],
            "n_infected_nbrs":  scores["n_infected_nbrs"],
            "degree":           deg,
            "betweenness":      round(b, 5),
            "composite_score":  round(composite, 5),
        })

    df = pd.DataFrame(rows).sort_values("composite_score", ascending=False).head(top_n)
    df = df.reset_index(drop=True)
    df.index += 1
    print(f"\n  Priority list (top {top_n}) — infected seeds: {infected}")
    print(df.to_string())
    return df


def plot_priority_list(df: pd.DataFrame, label: str):
    """Bar chart of top priority nodes."""
    fig, ax = plt.subplots(figsize=(9, 4))
    top = df.head(15)
    ax.barh(top["node"].astype(str), top["composite_score"],
            color="tomato", edgecolor="black", linewidth=0.4)
    ax.invert_yaxis()
    ax.set_xlabel("Composite Priority Score")
    ax.set_title(f"Priority List — Editors to Check Next\n({label})", fontsize=11)
    ax.grid(axis="x", alpha=0.3)
    plt.tight_layout()
    path = f"{OUTPUT_DIR}/task_c_priority_{label.lower().replace(' ','_')}.png"
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  [saved] {path}")


if __name__ == "__main__":
    from task_a import build_network, DATASETS

    for size, (name, path) in DATASETS.items():
        G, _, _ = build_network(path, f"{name} ({size})")
        Gc    = max(nx.connected_components(G), key=len)
        Gc    = G.subgraph(Gc).copy()
        nodes = list(Gc.nodes())

        print(f"\n{'='*55}")
        print(f"  Task C: {name} ({size})")
        print(f"{'='*55}")

        editor_a, editor_b = random.sample(nodes, 2)

        print("\n[Q1] Plausibility — behaviour NOT yet propagated:")
        propagation_plausibility(Gc, editor_a, editor_b)

        print("\n[Q2a] Priority list — ONE editor trolling:")
        df_one = priority_list(Gc, [editor_a], top_n=15)
        plot_priority_list(df_one, f"{name}_{size}_one_infected")

        print("\n[Q2b] Priority list — BOTH editors trolling:")
        df_both = priority_list(Gc, [editor_a, editor_b], top_n=15)
        plot_priority_list(df_both, f"{name}_{size}_both_infected")

    print("\nTask C complete.")
