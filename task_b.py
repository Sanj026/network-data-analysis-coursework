"""
Part 1 - Task B: Network Metrics
Situates each network in the regular / small-world / random continuum.
"""

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import kstest
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs_part1")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def largest_cc(G: nx.Graph) -> nx.Graph:
    """Return the largest connected component as a graph."""
    lcc = max(nx.connected_components(G), key=len)
    return G.subgraph(lcc).copy()


def erdos_renyi_reference(n: int, m: int, seed: int = 42) -> nx.Graph:
    """Erdős–Rényi random graph with same n, m as real network."""
    return nx.gnm_random_graph(n, m, seed=seed)


def compute_metrics(G: nx.Graph, label: str) -> dict:
    """Compute all Task B metrics for a graph."""
    Gc = largest_cc(G)
    n, m = Gc.number_of_nodes(), Gc.number_of_edges()
    degrees = [d for _, d in Gc.degree()]
    
    print(f"\n{'='*55}")
    print(f"  Metrics: {label}")
    print(f"{'='*55}")

    avg_degree  = np.mean(degrees)
    density     = nx.density(Gc)
    print(f"  LCC nodes / edges  : {n:,} / {m:,}")
    print(f"  Average degree     : {avg_degree:.4f}")
    print(f"  Density            : {density:.6f}")

    avg_clustering = nx.average_clustering(Gc)
    print(f"  Avg clustering coef: {avg_clustering:.4f}")

    if n > 2000:
        sample = list(Gc.nodes)[:500]
        sample_subG = Gc.subgraph(sample)
        sample_lcc = max(nx.connected_components(sample_subG), key=len)
        avg_path = nx.average_shortest_path_length(sample_subG.subgraph(sample_lcc))
        print(f"  Avg shortest path  : {avg_path:.4f}  (sampled 500 nodes)")
    else:
        avg_path = nx.average_shortest_path_length(Gc)
        print(f"  Avg shortest path  : {avg_path:.4f}")

    diameter = nx.diameter(Gc) if n <= 5000 else "skipped (too large)"
    print(f"  Diameter           : {diameter}")

    d_arr  = np.array(degrees)
    d_pos  = d_arr[d_arr > 0]
    ks_stat, ks_p = kstest(d_pos / d_pos.max(), 'powerlaw', args=(1,))
    print(f"  Degree KS stat (power-law): {ks_stat:.4f}  p={ks_p:.4f}")

    Gr = erdos_renyi_reference(n, m)
    Grc = largest_cc(Gr)
    c_rand = nx.average_clustering(Grc)
    if Grc.number_of_nodes() > 2000:
        rand_sample = Grc.subgraph(list(Grc.nodes)[:500])
        rand_lcc = max(nx.connected_components(rand_sample), key=len)
        l_rand = nx.average_shortest_path_length(rand_sample.subgraph(rand_lcc))
    else:
        l_rand = nx.average_shortest_path_length(Grc) if nx.is_connected(Grc) else float('nan')

    sigma = (avg_clustering / c_rand) / (avg_path / l_rand) if l_rand else float('nan')
    print(f"  ER random: C={c_rand:.4f}, L={l_rand:.4f}")
    print(f"  Small-world sigma  : {sigma:.4f}  (σ>1 → small-world)")

    if n <= 3000:
        bc = nx.betweenness_centrality(Gc, normalized=True)
    else:
        bc = nx.betweenness_centrality(Gc, normalized=True, k=300)
    top5_bc = sorted(bc.items(), key=lambda x: -x[1])[:5]
    print(f"  Top-5 betweenness centrality nodes: {[f'node {n} ({v:.4f})' for n,v in top5_bc]}")

    return {
        "label": label,
        "n": n, "m": m,
        "avg_degree": avg_degree,
        "density": density,
        "avg_clustering": avg_clustering,
        "avg_path": avg_path,
        "diameter": diameter,
        "sigma": sigma,
        "c_rand": c_rand,
        "l_rand": l_rand,
        "degrees": degrees,
        "betweenness": bc,
        "ks_stat": ks_stat,
        "ks_p": ks_p,
    }


def plot_metrics_summary(all_metrics: list):
    """One summary figure comparing all 3 networks."""
    fig = plt.figure(figsize=(16, 10))
    fig.suptitle("Network Metrics Summary — Wikidata Editor Networks",
                 fontsize=14, fontweight="bold")
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)

    colors = ["#3A86FF", "#FF6B6B", "#06D6A0"]

    for i, m in enumerate(all_metrics):
        ax = fig.add_subplot(gs[0, i])
        degrees = np.array(m["degrees"])
        degrees = degrees[degrees > 0]
        bins = np.logspace(np.log10(1), np.log10(degrees.max()+1), 40)
        ax.hist(degrees, bins=bins, color=colors[i], edgecolor="white",
                linewidth=0.3, density=True, alpha=0.85)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_title(f"{m['label']}\nDegree dist. (log-log)", fontsize=9)
        ax.set_xlabel("Degree"); ax.set_ylabel("Density")
        ax.grid(alpha=0.3)

    labels = [m["label"] for m in all_metrics]

    ax4 = fig.add_subplot(gs[1, 0])
    vals = [m["avg_clustering"] for m in all_metrics]
    rand = [m["c_rand"] for m in all_metrics]
    x = np.arange(len(labels))
    ax4.bar(x - 0.2, vals, 0.35, label="Real", color=colors)
    ax4.bar(x + 0.2, rand, 0.35, label="ER random", color="lightgrey", edgecolor="black")
    ax4.set_title("Avg Clustering Coefficient", fontsize=9)
    ax4.set_xticks(x); ax4.set_xticklabels(labels, fontsize=7)
    ax4.legend(fontsize=7); ax4.grid(axis="y", alpha=0.3)

    ax5 = fig.add_subplot(gs[1, 1])
    paths = [m["avg_path"] if isinstance(m["avg_path"], float) else 0 for m in all_metrics]
    l_rand = [m["l_rand"] if isinstance(m["l_rand"], float) else 0 for m in all_metrics]
    ax5.bar(x - 0.2, paths, 0.35, label="Real", color=colors)
    ax5.bar(x + 0.2, l_rand, 0.35, label="ER random", color="lightgrey", edgecolor="black")
    ax5.set_title("Avg Shortest Path Length", fontsize=9)
    ax5.set_xticks(x); ax5.set_xticklabels(labels, fontsize=7)
    ax5.legend(fontsize=7); ax5.grid(axis="y", alpha=0.3)

    ax6 = fig.add_subplot(gs[1, 2])
    sigmas = [m["sigma"] for m in all_metrics]
    bars = ax6.bar(labels, sigmas, color=colors, edgecolor="black", linewidth=0.5)
    ax6.axhline(1, color="red", linestyle="--", linewidth=1, label="σ=1 (threshold)")
    ax6.set_title("Small-World Sigma (σ)", fontsize=9)
    ax6.set_xticklabels(labels, fontsize=7)
    ax6.legend(fontsize=7); ax6.grid(axis="y", alpha=0.3)
    for bar, v in zip(bars, sigmas):
        ax6.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                 f"{v:.2f}", ha="center", fontsize=8)

    plt.savefig(f"{OUTPUT_DIR}/task_b_metrics.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  [saved] {OUTPUT_DIR}/task_b_metrics.png")


if __name__ == "__main__":
    from task_a import build_network, DATASETS

    all_metrics = []
    graphs = {}
    for size, (name, path) in DATASETS.items():
        G, _, _ = build_network(path, f"{name} ({size})")
        short_label = size.capitalize()
        m = compute_metrics(G, short_label)
        all_metrics.append(m)

    plot_metrics_summary(all_metrics)
    print("\nTask B complete.")
