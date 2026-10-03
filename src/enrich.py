import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gseapy

COLUMNS = ["term", "library", "direction", "adj_p", "overlap", "genes"]
MIN_DIRECTION = 20


def direction_lists(table, set_name):
    sub = table[table["set"] == set_name]
    fc = sub["aging_logFC"] if set_name == "aging_only" else sub["ad_logFC"]
    up = sub.loc[fc > 0, "gene"].tolist()
    down = sub.loc[fc < 0, "gene"].tolist()
    if len(up) >= MIN_DIRECTION and len(down) >= MIN_DIRECTION:
        return {"up": up, "down": down}
    return {"all": sub["gene"].tolist()}


def enrichr_one(genes, gene_sets, direction):
    res = gseapy.enrichr(gene_list=genes, gene_sets=gene_sets, organism="human", outdir=None, no_plot=True).results
    return pd.DataFrame({
        "term": res["Term"],
        "library": res["Gene_set"],
        "direction": direction,
        "adj_p": res["Adjusted P-value"],
        "overlap": res["Overlap"],
        "genes": res["Genes"],
    })


def enrich_set(table, set_name, gene_sets):
    frames = [enrichr_one(genes, gene_sets, direction) for direction, genes in direction_lists(table, set_name).items() if genes]
    if not frames:
        return pd.DataFrame(columns=COLUMNS)
    return pd.concat(frames).sort_values("adj_p").reset_index(drop=True)


def plot_enrichment(result, set_name, path):
    top = result.head(10).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 5))
    if top.empty:
        ax.text(0.5, 0.5, "No enriched terms", ha="center", va="center")
        ax.axis("off")
    else:
        labels = [f"{t[:60]} [{d}]" for t, d in zip(top["term"], top["direction"])]
        colors = top["direction"].map({"up": "#c0392b", "down": "#2c7bb6", "all": "#7f8c8d"})
        ax.barh(labels, -np.log10(top["adj_p"]), color=colors)
        ax.set_xlabel("-log10 adjusted p")
    ax.set_title(f"Top pathways: {set_name}")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_enrich(table, set_name, gene_sets, results_dir, figures_dir, skip=False):
    path = results_dir / f"enrichment_{set_name}.csv"
    result = pd.DataFrame(columns=COLUMNS)
    try:
        if not skip:
            result = enrich_set(table, set_name, gene_sets)
    except Exception as err:
        print(f"  Enrichr unreachable for {set_name}: {type(err).__name__}; writing empty table")
    result.to_csv(path, index=False)
    plot_enrichment(result, set_name, figures_dir / f"enrichment_{set_name}.png")
    return result
