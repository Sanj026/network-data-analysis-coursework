"""
Part 1 - Task A: Network Construction
Wikidata Editor Social Network
Files: BOT_REQUESTS (small), PROJECT_CHAT (medium), REQUEST_FOR_DELETION (large)
"""

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations
from collections import defaultdict
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

DATASETS = {
    "small":  ("BOT_REQUESTS",         os.path.join(BASE_DIR, "datasets", "BOT_REQUESTS.csv")),
    "medium": ("PROJECT_CHAT",         os.path.join(BASE_DIR, "datasets", "PROJECT_CHAT.csv")),
    "large":  ("REQUEST_FOR_DELETION", os.path.join(BASE_DIR, "datasets", "REQUEST_FOR_DELETION.csv")),
}
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs_part1")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Core builder ──────────────────────────────────────────────────────────────
def build_network(filepath: str, label: str) -> tuple[nx.Graph, dict, pd.DataFrame]:
    """
    Build a co-comment editor network from a CSV file.

    Nodes  : unique users (integer node ID)
    Edges  : pair of users who both commented in the SAME (page_name, thread_subject)
    Returns: (Graph, node_metadata dict, raw dataframe)
    """
    df = pd.read_csv(filepath)
    df.columns = df.columns.str.strip()
    df = df.dropna(subset=["username", "page_name", "thread_subject"])
    df["username"]       = df["username"].str.strip()
    df["page_name"]      = df["page_name"].str.strip()
    df["thread_subject"] = df["thread_subject"].str.strip()

    unique_users = sorted(df["username"].unique())
    user_to_id   = {u: i for i, u in enumerate(unique_users)}

    node_metadata = {}
    for user, grp in df.groupby("username"):
        node_metadata[user_to_id[user]] = {
            "username":      user,
            "n_comments":    len(grp),
            "n_pages":       grp["page_name"].nunique(),
            "n_threads":     grp["thread_subject"].nunique(),
        }

    G = nx.Graph()
    G.add_nodes_from(range(len(unique_users)))
    nx.set_node_attributes(G, node_metadata)

    thread_groups = (
        df.groupby(["page_name", "thread_subject"])["username"]
        .apply(set)
        .reset_index()
    )

    edge_weight = defaultdict(int)
    for _, row in thread_groups.iterrows():
        users_in_thread = list(row["username"])
        for u, v in combinations(users_in_thread, 2):
            uid, vid = user_to_id[u], user_to_id[v]
            key = (min(uid, vid), max(uid, vid))
            edge_weight[key] += 1

    for (u, v), w in edge_weight.items():
        G.add_edge(u, v, weight=w)

    print(f"\n{'='*55}")
    print(f"  Network: {label}")
    print(f"{'='*55}")
    print(f"  Rows in CSV        : {len(df):,}")
    print(f"  Unique users (nodes): {G.number_of_nodes():,}")
    print(f"  Edges              : {G.number_of_edges():,}")
    print(f"  Threads (groups)   : {len(thread_groups):,}")
    print(f"  Connected components: {nx.number_connected_components(G):,}")

    return G, node_metadata, df


def plot_network(G: nx.Graph, label: str, max_nodes: int = 300):
    """Visualise a (sampled) subgraph with node size ∝ degree."""
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_title(f"Wikidata Editor Network — {label}\n"
                 f"({G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges)",
                 fontsize=13, fontweight="bold")

    lcc = max(nx.connected_components(G), key=len)
    subG = G.subgraph(lcc).copy()
    if subG.number_of_nodes() > max_nodes:
        sampled = list(subG.nodes)[:max_nodes]
        subG = subG.subgraph(sampled).copy()

    degrees = dict(subG.degree())
    node_sizes = [max(20, degrees[n] * 15) for n in subG.nodes()]
    node_colors = [degrees[n] for n in subG.nodes()]

    pos = nx.spring_layout(subG, seed=42, k=0.5)
    nx.draw_networkx_edges(subG, pos, ax=ax, alpha=0.2, width=0.5, edge_color="grey")
    sc = nx.draw_networkx_nodes(subG, pos, ax=ax,
                                 node_size=node_sizes,
                                 node_color=node_colors,
                                 cmap=plt.cm.plasma, alpha=0.85)
    plt.colorbar(sc, ax=ax, label="Degree")
    ax.axis("off")
    plt.tight_layout()
    safe_label = label.lower().replace(" ", "_").replace("(", "").replace(")", "")
    path = f"{OUTPUT_DIR}/network_{safe_label}.png"
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  [saved] {path}")


def plot_degree_distribution(graphs: dict):
    """Plot degree distributions for all 3 networks side by side."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    fig.suptitle("Degree Distributions — Wikidata Editor Networks", fontsize=13, fontweight="bold")

    for ax, (size, (name, G)) in zip(axes, graphs.items()):
        degrees = [d for _, d in G.degree()]
        ax.hist(degrees, bins=50, color="steelblue", edgecolor="white", linewidth=0.4, log=True)
        ax.set_title(f"{name}\n({size})", fontsize=10)
        ax.set_xlabel("Degree")
        ax.set_ylabel("Count (log scale)")
        ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    path = f"{OUTPUT_DIR}/degree_distributions.png"
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"  [saved] {path}")


if __name__ == "__main__":
    graphs = {}

    for size, (name, path) in DATASETS.items():
        G, _, _ = build_network(path, f"{name} ({size})")
        graphs[size] = (name, G)
        plot_network(G, f"{name} ({size})")

    plot_degree_distribution(graphs)
    print("\nTask A complete. Outputs saved to outputs/")
