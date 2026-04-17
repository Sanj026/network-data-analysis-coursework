"""
Part 2 - Leeds Road Network Analysis (Tasks A, B, C)

Runs end-to-end:
1) Builds a 1 km^2 drivable network around Leeds city centre area with high accident density.
2) Computes spatial network metrics, circuitry, and planarity checks.
3) Analyses accident spatial patterns, Moran's I, K-function, and intersection proximity.
4) Builds network Voronoi cells (N=4) and searches approximate 42 km closed loops.

Outputs:
- PNG figures in outputs_part2/
- JSON summary in outputs_part2/part2_summary.json
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict, List, Tuple

# Ensure matplotlib cache is writable in workspace before importing pyplot
os.environ.setdefault("MPLCONFIGDIR", str(Path(".mplconfig").resolve()))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import osmnx as ox
import pandas as pd
import geopandas as gpd
import spaghetti
from esda.moran import Moran
from libpysal.weights import W
from pointpats import k_test
from pyproj import Transformer


BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = (BASE_DIR / "datasets").resolve()
OUTPUT_DIR = (BASE_DIR / "outputs_part2").resolve()
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

LEEDS_CENTRE_LAT = 53.7997
LEEDS_CENTRE_LON = -1.5492
SQUARE_SIDE_M = 1000
YEARS_WITH_COORDS = list(range(2009, 2017))


def load_accidents() -> pd.DataFrame:
    """Load and harmonise accident files that include Easting/Northing."""
    frames = []
    for year in YEARS_WITH_COORDS:
        path = DATA_DIR / f"{year}.csv"
        if not path.exists():
            continue
        df = pd.read_csv(path, encoding="latin1", low_memory=False)
        cols = {c.strip(): c for c in df.columns}
        east_col = cols.get("Grid Ref: Easting", cols.get("Easting"))
        north_col = cols.get("Grid Ref: Northing", cols.get("Northing"))
        if east_col is None or north_col is None:
            continue
        sub = df[[east_col, north_col]].copy()
        sub.columns = ["easting", "northing"]
        sub["year"] = year
        frames.append(sub)

    if not frames:
        raise FileNotFoundError("No accident files with coordinates were found.")

    out = pd.concat(frames, ignore_index=True)
    out["easting"] = pd.to_numeric(out["easting"], errors="coerce")
    out["northing"] = pd.to_numeric(out["northing"], errors="coerce")
    out = out.dropna(subset=["easting", "northing"]).copy()
    out = out[(out["easting"] > 0) & (out["northing"] > 0)]
    out["accident_id"] = np.arange(len(out))
    return out


def choose_1km_square(acc: pd.DataFrame) -> Tuple[Tuple[float, float], Tuple[float, float, float, float], int]:
    """
    Search around Leeds centre for a 1 km^2 square containing many accidents.
    Returns (centre_bng, bbox_bng, accident_count).
    """
    to_bng = Transformer.from_crs(4326, 27700, always_xy=True)
    cx, cy = to_bng.transform(LEEDS_CENTRE_LON, LEEDS_CENTRE_LAT)
    half = SQUARE_SIDE_M / 2

    best_count = -1
    best_center = (cx, cy)
    for dx in range(-1500, 1501, 250):
        for dy in range(-1500, 1501, 250):
            x = cx + dx
            y = cy + dy
            xmin, xmax = x - half, x + half
            ymin, ymax = y - half, y + half
            mask = (
                (acc["easting"] >= xmin)
                & (acc["easting"] <= xmax)
                & (acc["northing"] >= ymin)
                & (acc["northing"] <= ymax)
            )
            count = int(mask.sum())
            if count > best_count:
                best_count = count
                best_center = (x, y)

    bx, by = best_center
    bbox = (bx - half, by - half, bx + half, by + half)
    return best_center, bbox, best_count


def bbox_to_latlon(bbox_bng: Tuple[float, float, float, float]) -> Tuple[float, float, float, float]:
    """Convert BNG bbox (xmin,ymin,xmax,ymax) to (north,south,east,west)."""
    xmin, ymin, xmax, ymax = bbox_bng
    to_wgs = Transformer.from_crs(27700, 4326, always_xy=True)
    west, south = to_wgs.transform(xmin, ymin)
    east, north = to_wgs.transform(xmax, ymax)
    return north, south, east, west


def build_network(north: float, south: float, east: float, west: float) -> nx.MultiDiGraph:
    return ox.graph_from_bbox((west, south, east, north), network_type="drive", simplify=True)


def edge_key(u: int, v: int, k: int) -> str:
    return f"{u}|{v}|{k}"


def compute_edge_accident_mapping(G_proj: nx.MultiDiGraph, acc_area: pd.DataFrame) -> Dict[str, int]:
    """
    Snap accidents to nearest edges and count accidents per directed edge key.
    """
    xs, ys = acc_area["easting"].to_numpy(), acc_area["northing"].to_numpy()
    nearest = ox.distance.nearest_edges(G_proj, X=xs, Y=ys)

    counts: Dict[str, int] = {}
    for u, v, k in nearest:
        ek = edge_key(u, v, k)
        counts[ek] = counts.get(ek, 0) + 1
    return counts


def moran_on_edges(G_u: nx.Graph, edge_counts: Dict[str, int]) -> Tuple[float, float]:
    """
    Build edge adjacency graph and compute Moran's I on accident counts.
    """
    undirected_edges = list(G_u.edges())
    idx = {e: i for i, e in enumerate(undirected_edges)}
    y = np.zeros(len(undirected_edges), dtype=float)

    for (u, v) in undirected_edges:
        directed_sum = 0
        for k in G_u[u][v]:
            directed_sum += edge_counts.get(edge_key(u, v, k), 0)
            directed_sum += edge_counts.get(edge_key(v, u, k), 0)
        y[idx[(u, v)]] = directed_sum

    neighbors = {i: set() for i in range(len(undirected_edges))}
    node_to_edges: Dict[int, List[int]] = {}
    for i, (u, v) in enumerate(undirected_edges):
        node_to_edges.setdefault(u, []).append(i)
        node_to_edges.setdefault(v, []).append(i)

    for eids in node_to_edges.values():
        for i in eids:
            for j in eids:
                if i != j:
                    neighbors[i].add(j)
    neighbors = {k: list(v) for k, v in neighbors.items()}
    w = W(neighbors)
    mi = Moran(y, w, permutations=999)
    return float(mi.I), float(mi.p_sim)


def estimate_planarity_crossings(G_u: nx.Graph) -> int:
    """
    Approximate non-planarity by counting geometric line crossings without shared nodes.
    Uses a subset for tractability.
    """
    gdf_edges = ox.graph_to_gdfs(G_u, nodes=False, edges=True).reset_index()
    sample = gdf_edges.sample(n=min(len(gdf_edges), 400), random_state=42)
    geoms = list(sample["geometry"])
    uv = list(zip(sample["u"], sample["v"]))
    crossings = 0
    for i in range(len(geoms)):
        for j in range(i + 1, len(geoms)):
            if set(uv[i]) & set(uv[j]):
                continue
            if geoms[i].crosses(geoms[j]):
                crossings += 1
    return crossings


def nearest_intersection_fraction(
    G_u: nx.Graph, acc_area: pd.DataFrame
) -> Tuple[np.ndarray, float]:
    edges = ox.graph_to_gdfs(G_u, nodes=False, edges=True)[["geometry"]].reset_index(drop=True)
    ntw = spaghetti.Network(in_data=edges)
    points = gpd.GeoDataFrame(
        acc_area[["accident_id"]].copy(),
        geometry=gpd.points_from_xy(acc_area["easting"], acc_area["northing"]),
        crs="EPSG:27700",
    )
    ntw.snapobservations(points, "accidents", idvariable="accident_id", attribute=True)
    snapped = ntw.pointpatterns["accidents"]

    fracs = []
    for dist_map in snapped.dist_to_vertex.values():
        dists = list(dist_map.values())
        total = sum(dists)
        fracs.append(min(dists) / total if total else 0.0)
    arr = np.array(fracs, dtype=float)
    return arr, float(arr.mean())


def build_closed_route(G_u: nx.Graph, seed: int, waypoints: np.ndarray) -> Tuple[List[int], float] | None:
    segments = []
    current = seed
    for waypoint in waypoints:
        try:
            seg = nx.shortest_path(G_u, current, int(waypoint), weight="length")
        except nx.NetworkXNoPath:
            return None
        segments.append(seg)
        current = int(waypoint)
    try:
        segments.append(nx.shortest_path(G_u, current, seed, weight="length"))
    except nx.NetworkXNoPath:
        return None

    route = segments[0]
    for seg in segments[1:]:
        route.extend(seg[1:])
    return route, float(nx.path_weight(G_u, route, weight="length"))


def repeat_route(route: List[int], length: float, target_m: float, tol_m: float) -> Tuple[List[int], float] | None:
    for mult in range(1, 5):
        total = length * mult
        if target_m - tol_m <= total <= target_m + tol_m:
            repeated = route[:]
            for _ in range(mult - 1):
                repeated.extend(route[1:])
            return repeated, total
    return None


def search_cell_loops(
    G_u: nx.Graph,
    seed: int,
    cell_nodes: List[int],
    target_m: float = 42000,
    tol_m: float = 5000,
    tries: int = 1200,
) -> List[Tuple[List[int], float]]:
    if len(cell_nodes) < 12:
        return []

    rng = np.random.default_rng(42 + int(seed) % 1000)
    pool = np.array(cell_nodes)
    waypoint_counts = (4, 5, 6, 7, 8)
    candidates: List[Tuple[List[int], float]] = []

    for count in waypoint_counts:
        if len(pool) <= count:
            continue
        for _ in range(tries // len(waypoint_counts)):
            waypoints = rng.choice(pool, size=count, replace=False)
            built = build_closed_route(G_u, seed, waypoints)
            if built is None:
                continue
            route, length = built
            repeated = repeat_route(route, length, target_m, tol_m)
            if repeated is not None:
                candidates.append(repeated)

    dedup = {}
    for route, length in candidates:
        key = tuple(route[:20]) + (len(route), round(length, 1))
        dedup[key] = (route, float(length))
    return sorted(dedup.values(), key=lambda item: abs(item[1] - target_m))[:3]


def shifted_seed_candidates(G_u: nx.Graph, seed: int, cell_nodes: List[int], max_candidates: int = 3) -> List[int]:
    node_gdf = ox.graph_to_gdfs(G_u, nodes=True, edges=False)
    sx, sy = float(node_gdf.loc[seed, "x"]), float(node_gdf.loc[seed, "y"])
    ranked = []
    for node in cell_nodes:
        if node == seed:
            continue
        x, y = float(node_gdf.loc[node, "x"]), float(node_gdf.loc[node, "y"])
        ranked.append((((x - sx) ** 2 + (y - sy) ** 2), -G_u.degree(node), int(node)))
    ranked.sort()
    return [node for _, _, node in ranked[:max_candidates]]


def subdivide_cell(G_u: nx.Graph, cell_nodes: List[int]) -> List[List[int]]:
    node_gdf = ox.graph_to_gdfs(G_u, nodes=True, edges=False).loc[cell_nodes]
    x_mid = float(node_gdf["x"].median())
    y_mid = float(node_gdf["y"].median())
    groups = [
        node_gdf[node_gdf["x"] <= x_mid].index.tolist(),
        node_gdf[node_gdf["x"] > x_mid].index.tolist(),
        node_gdf[node_gdf["y"] <= y_mid].index.tolist(),
        node_gdf[node_gdf["y"] > y_mid].index.tolist(),
    ]
    return [group for group in groups if len(group) >= 12]


def find_42km_cycles(
    G_u: nx.Graph,
    seed_to_nodes: Dict[int, List[int]],
    target_m: float = 42000,
    tol_m: float = 5000,
    tries_per_seed: int = 1200,
) -> Tuple[Dict[int, List[Tuple[List[int], float]]], Dict[int, str]]:
    results: Dict[int, List[Tuple[List[int], float]]] = {}
    notes: Dict[int, str] = {}

    for seed, cell_nodes in seed_to_nodes.items():
        attempts = ["base search"]
        candidates = search_cell_loops(G_u, seed, cell_nodes, target_m=target_m, tol_m=tol_m, tries=tries_per_seed)

        if not candidates:
            attempts.append("shifted seed")
            for alt_seed in shifted_seed_candidates(G_u, seed, cell_nodes):
                candidates = search_cell_loops(G_u, alt_seed, cell_nodes, target_m=target_m, tol_m=tol_m, tries=tries_per_seed)
                if candidates:
                    break

        if not candidates:
            attempts.append("subdivided cell")
            for group in subdivide_cell(G_u, cell_nodes):
                subgroup_seed = seed if seed in group else group[0]
                candidates = search_cell_loops(
                    G_u,
                    subgroup_seed,
                    group,
                    target_m=target_m,
                    tol_m=tol_m,
                    tries=max(400, tries_per_seed // 2),
                )
                if candidates:
                    break

        results[seed] = candidates
        notes[seed] = f"{'succeeded' if candidates else 'failed'}; attempts: {', '.join(attempts)}"

    return results, notes


def build_seed_cells(G_u: nx.Graph, center_xy: Tuple[float, float]) -> Tuple[List[int], Dict[int, List[int]]]:
    """Create 4 seeds (quadrants, high-degree) and assign each node to nearest seed."""
    nodes_gdf = ox.graph_to_gdfs(G_u, nodes=True, edges=False)
    cx, cy = center_xy
    candidates = nodes_gdf.copy()
    candidates["deg"] = [G_u.degree(n) for n in candidates.index]
    quadrants = [
        (candidates["x"] <= cx) & (candidates["y"] <= cy),
        (candidates["x"] > cx) & (candidates["y"] <= cy),
        (candidates["x"] <= cx) & (candidates["y"] > cy),
        (candidates["x"] > cx) & (candidates["y"] > cy),
    ]
    seed_nodes = []
    for q in quadrants:
        qdf = candidates[q]
        if qdf.empty:
            qdf = candidates
        seed_nodes.append(int(qdf.sort_values("deg", ascending=False).index[0]))

    cells = {s: [] for s in seed_nodes}
    for n in G_u.nodes():
        dists = []
        for s in seed_nodes:
            try:
                d = nx.shortest_path_length(G_u, s, n, weight="length")
            except nx.NetworkXNoPath:
                d = np.inf
            dists.append((s, d))
        best = min(dists, key=lambda x: x[1])[0]
        cells[best].append(n)
    return seed_nodes, cells


def plot_marathon_loops(
    G_u: nx.Graph,
    loop_results: Dict[int, List[Tuple[List[int], float]]],
    output_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 8))
    ox.plot_graph(G_u, ax=ax, node_size=0, edge_color="lightgrey", edge_linewidth=0.5, show=False, close=False)
    colors = ["tomato", "royalblue", "seagreen", "darkorange"]
    drawn = 0

    for color, (_, routes) in zip(colors, loop_results.items()):
        if not routes:
            continue
        route = routes[0][0]
        xs, ys = [], []
        for node in route:
            xs.append(G_u.nodes[node]["x"])
            ys.append(G_u.nodes[node]["y"])
        ax.plot(xs, ys, color=color, linewidth=2, alpha=0.9)
        drawn += 1

    ax.set_title(f"Marathon loops found across cells ({drawn} cells)")
    plt.tight_layout()
    plt.savefig(output_path, dpi=160)
    plt.close()


def plot_voronoi_map(G_u: nx.Graph, cells: Dict[int, List[int]], output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 8))
    ox.plot_graph(G_u, ax=ax, node_size=0, edge_color="lightgrey", edge_linewidth=0.4, show=False, close=False)
    colors = ["tomato", "royalblue", "seagreen", "darkorange"]
    for color, (seed, nodes) in zip(colors, cells.items()):
        xs = [G_u.nodes[node]["x"] for node in nodes]
        ys = [G_u.nodes[node]["y"] for node in nodes]
        ax.scatter(xs, ys, s=6, color=color, alpha=0.6)
        ax.scatter(G_u.nodes[seed]["x"], G_u.nodes[seed]["y"], s=60, color=color, edgecolor="black")
    ax.set_title("Network Voronoi partition (N=4)")
    plt.tight_layout()
    plt.savefig(output_path, dpi=160)
    plt.close()


def main() -> None:
    print("Loading accidents...")
    accidents = load_accidents()

    print("Selecting best 1 km^2 area near Leeds centre...")
    center_bng, bbox_bng, n_acc = choose_1km_square(accidents)
    xmin, ymin, xmax, ymax = bbox_bng
    north, south, east, west = bbox_to_latlon(bbox_bng)

    acc_area = accidents[
        (accidents["easting"] >= xmin)
        & (accidents["easting"] <= xmax)
        & (accidents["northing"] >= ymin)
        & (accidents["northing"] <= ymax)
    ].copy()

    print(f"Area accident count: {len(acc_area)}")
    print("Downloading OSM road graph...")
    G = build_network(north, south, east, west)
    G_proj = ox.project_graph(G, to_crs="EPSG:27700")
    G_u = ox.convert.to_undirected(G_proj)

    area_km2 = 1.0
    stats = ox.stats.basic_stats(G_proj, area=area_km2 * 1_000_000)
    diameter = nx.diameter(max((G_u.subgraph(c) for c in nx.connected_components(G_u)), key=len))
    avg_street_length = float(stats.get("street_length_avg", np.nan))
    node_density = float(stats.get("node_density_km", np.nan))
    intersection_density = float(stats.get("intersection_density_km", np.nan))
    edge_density = float(stats.get("edge_density_km", np.nan))
    circuity = float(stats.get("circuity_avg", np.nan))
    crossings_est = estimate_planarity_crossings(G_u)

    edge_counts = compute_edge_accident_mapping(G_proj, acc_area)
    moran_i, moran_p = moran_on_edges(G_u, edge_counts)

    coords = acc_area[["easting", "northing"]].to_numpy()
    support = np.linspace(50, 500, 10)
    k = k_test(coords, support=support, keep_simulations=False, n_simulations=99, n_jobs=1)
    k_pvalue_min = float(np.min(k.pvalue))

    frac_arr, frac_mean = nearest_intersection_fraction(G_u, acc_area)

    seed_nodes, cells = build_seed_cells(G_u, center_bng)
    cycle_results, cycle_notes = find_42km_cycles(G_u, cells, target_m=42000, tol_m=5000, tries_per_seed=1200)
    feasible_cells = sum(1 for s in seed_nodes if len(cycle_results[s]) > 0)

    ref_half = 4000
    ref_bbox = (center_bng[0] - ref_half, center_bng[1] - ref_half, center_bng[0] + ref_half, center_bng[1] + ref_half)
    ref_n, ref_s, ref_e, ref_w = bbox_to_latlon(ref_bbox)
    G_ref = build_network(ref_n, ref_s, ref_e, ref_w)
    G_ref_proj = ox.project_graph(G_ref, to_crs="EPSG:27700")
    G_ref_u = ox.convert.to_undirected(G_ref_proj)
    ref_seeds, ref_cells = build_seed_cells(G_ref_u, center_bng)
    ref_cycles, ref_notes = find_42km_cycles(G_ref_u, ref_cells, target_m=42000, tol_m=5000, tries_per_seed=1400)
    ref_feasible = sum(1 for s in ref_seeds if len(ref_cycles[s]) > 0)

    plt.figure(figsize=(8, 8))
    ox.plot_graph(G, node_size=3, edge_linewidth=0.6, show=False, close=False)
    plt.title("Leeds 1 km² road network")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "part2_road_network.png", dpi=160)
    plt.close()

    fig, ax = plt.subplots(figsize=(8, 8))
    ox.plot_graph(G, ax=ax, node_size=3, edge_linewidth=0.6, show=False, close=False)
    to_wgs = Transformer.from_crs(27700, 4326, always_xy=True)
    lon, lat = to_wgs.transform(acc_area["easting"].to_numpy(), acc_area["northing"].to_numpy())
    ax.scatter(lon, lat, s=8, c="red", alpha=0.5, label="Accidents")
    ax.legend()
    ax.set_title("Leeds 1 km² area road network + accidents")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "part2_area_accidents.png", dpi=160)
    plt.close()

    plt.figure(figsize=(7, 4))
    plt.hist(frac_arr, bins=25, color="teal", edgecolor="white")
    plt.xlabel("Fraction of edge length to nearest intersection (0=at junction)")
    plt.ylabel("Accident count")
    plt.title("Accident proximity to intersections")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "part2_intersection_fraction.png", dpi=160)
    plt.close()

    voronoi_graph = G_ref_u if ref_feasible >= feasible_cells else G_u
    voronoi_cells = ref_cells if ref_feasible >= feasible_cells else cells
    plot_voronoi_map(voronoi_graph, voronoi_cells, OUTPUT_DIR / "part2_voronoi_map.png")

    loop_plot_results = ref_cycles if ref_feasible >= feasible_cells else cycle_results
    loop_plot_graph = G_ref_u if ref_feasible >= feasible_cells else G_u
    plot_marathon_loops(loop_plot_graph, loop_plot_results, OUTPUT_DIR / "part2_marathon_loops.png")

    loop_results_note = {
        "initial": {
            "succeeded_cells": [int(s) for s in seed_nodes if cycle_results[s]],
            "failed_cells": [int(s) for s in seed_nodes if not cycle_results[s]],
            "notes": {str(s): cycle_notes[s] for s in seed_nodes},
        },
        "refinement": {
            "succeeded_cells": [int(s) for s in ref_seeds if ref_cycles[s]],
            "failed_cells": [int(s) for s in ref_seeds if not ref_cycles[s]],
            "notes": {str(s): ref_notes[s] for s in ref_seeds},
        },
    }

    summary = {
        "selected_area": {
            "bbox_bng": {"xmin": xmin, "ymin": ymin, "xmax": xmax, "ymax": ymax},
            "bbox_wgs84": {"north": north, "south": south, "east": east, "west": west},
            "centre_bng": {"x": center_bng[0], "y": center_bng[1]},
            "area_km2": area_km2,
            "accidents_in_area": int(len(acc_area)),
        },
        "task_a_metrics": {
            "diameter_nodes": int(diameter),
            "average_street_length_m": avg_street_length,
            "node_density_per_km2": node_density,
            "intersection_density_per_km2": intersection_density,
            "edge_density_m_per_km2": edge_density,
            "average_circuity": circuity,
            "estimated_nonplanar_crossings_sample": int(crossings_est),
        },
        "task_b_metrics": {
            "morans_i": moran_i,
            "morans_i_pvalue": moran_p,
            "k_function_min_pvalue": k_pvalue_min,
            "mean_intersection_fraction": frac_mean,
        },
        "task_c_metrics": {
            "initial_seed_nodes": seed_nodes,
            "initial_feasible_cells_for_42km_loop": int(feasible_cells),
            "initial_cycle_lengths_m": {
                str(s): [round(length, 2) for _, length in cycle_results[s]]
                for s in seed_nodes
            },
            "refinement": {
                "area": "8km x 8km around centre",
                "seed_nodes": ref_seeds,
                "feasible_cells_for_42km_loop": int(ref_feasible),
                "cycle_lengths_m": {
                    str(s): [round(length, 2) for _, length in ref_cycles[s]]
                    for s in ref_seeds
                },
            },
            "results_note": loop_results_note,
            "final_loop_count": int(max(feasible_cells, ref_feasible)),
        },
    }

    with open(OUTPUT_DIR / "part2_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("Saved outputs to", OUTPUT_DIR)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
