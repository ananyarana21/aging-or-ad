import numpy as np
import pandas as pd
import networkx as nx
import requests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SET_GROUPS = {
    "shared": ["young", "old", "AD"],
    "aging_only": ["young", "old"],
    "ad_only": ["old", "AD"],
    "ad": ["old", "AD"],
}
STRING_URL = "https://string-db.org/api/json/network"


def coexpr_graph(expr, meta, genes, set_name, r_cut):
    samples = meta.loc[meta["group"].isin(SET_GROUPS[set_name]), "sample"]
    sub = expr.loc[genes, samples]
    sub = sub[sub.std(axis=1) > 0]
    graph = nx.Graph()
    graph.add_nodes_from(sub.index)
    if len(sub) < 2:
        return graph
    corr = np.corrcoef(sub.to_numpy())
    rows, cols = np.where(np.triu(np.abs(corr) >= r_cut, k=1))
    graph.add_edges_from(zip(sub.index[rows], sub.index[cols]))
    return graph


def string_degrees(genes):
    try:
        resp = requests.post(STRING_URL, data={"identifiers": "\r".join(genes), "species": 9606, "caller_identity": "brain-aging-vs-ad"}, timeout=120)
        resp.raise_for_status()
        edges = pd.DataFrame(resp.json())
    except Exception as err:
        print(f"  STRING unreachable: {type(err).__name__}")
        return {}
    if edges.empty:
        return {}
    graph = nx.from_pandas_edgelist(edges, "preferredName_A", "preferredName_B")
    return dict(graph.degree())


def top_hubs(graph, table, set_name, top_n, use_string):
    degrees = sorted(graph.degree(), key=lambda x: (-x[1], x[0]))[:top_n]
    hubs = pd.DataFrame(degrees, columns=["gene", "degree"])
    hubs.insert(0, "set", set_name)
    fc = table.set_index("gene")[["aging_logFC", "ad_logFC"]]
    hubs["logFC_aging"] = hubs["gene"].map(fc["aging_logFC"])
    hubs["logFC_ad"] = hubs["gene"].map(fc["ad_logFC"])
    if use_string:
        string_deg = string_degrees(list(graph.nodes))
        hubs["string_degree"] = hubs["gene"].map(string_deg).fillna(0).astype(int)
    return hubs


def label_hubs(ax, pos, hub_set):
    coords = np.array(list(pos.values()))
    center = coords.mean(axis=0)
    radius = np.abs(coords - center).max() * 1.15
    hubs = sorted(hub_set, key=lambda n: np.arctan2(*(pos[n] - center)[::-1]))
    for i, name in enumerate(hubs):
        angle = 2 * np.pi * i / len(hubs)
        target = center + radius * np.array([np.cos(angle), np.sin(angle)])
        ax.annotate(name, xy=pos[name], xytext=target, fontsize=8, ha="center", va="center",
                    bbox={"fc": "white", "ec": "#c0392b", "lw": 0.5, "pad": 1.5},
                    arrowprops={"arrowstyle": "-", "color": "#c0392b", "lw": 0.5})
    pad = radius * 1.25
    ax.set_xlim(center[0] - pad, center[0] + pad)
    ax.set_ylim(center[1] - pad, center[1] + pad)


def plot_network(graph, hubs, set_name, path, seed):
    fig, ax = plt.subplots(figsize=(7, 7))
    if graph.number_of_edges() == 0:
        ax.text(0.5, 0.5, "No edges at this correlation cutoff", ha="center", va="center")
    else:
        lcc = graph.subgraph(max(nx.connected_components(graph), key=len))
        pos = nx.spring_layout(lcc, seed=seed, k=4 / np.sqrt(len(lcc)), iterations=200)
        deg = dict(lcc.degree())
        sizes = [10 + 200 * deg[n] / max(deg.values()) for n in lcc.nodes]
        hub_set = set(hubs["gene"]) & set(lcc.nodes)
        colors = ["#c0392b" if n in hub_set else "#7f8c8d" for n in lcc.nodes]
        nx.draw_networkx_edges(lcc, pos, ax=ax, alpha=0.15, width=0.4)
        nx.draw_networkx_nodes(lcc, pos, ax=ax, node_size=sizes, node_color=colors, linewidths=0)
        label_hubs(ax, pos, hub_set)
        ax.set_title(f"{set_name}: largest component ({lcc.number_of_nodes()} genes, {lcc.number_of_edges()} edges)")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_network(expr, meta, table, set_names, cfg, figures_dir, use_string):
    frames = []
    for set_name in set_names:
        genes = [g for g in table.loc[table["set"] == set_name, "gene"] if g in expr.index]
        graph = coexpr_graph(expr, meta, genes, set_name, cfg.COEXPR_R)
        hubs = top_hubs(graph, table, set_name, cfg.TOP_HUBS, use_string)
        plot_network(graph, hubs, set_name, figures_dir / f"network_{set_name}.png", cfg.SEED)
        frames.append(hubs)
    return pd.concat(frames, ignore_index=True)
