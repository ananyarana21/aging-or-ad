import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from scipy import stats

SET_COLORS = {"shared": "#8e44ad", "aging_only": "#2c7bb6", "ad_only": "#c0392b"}


def assign_set(row):
    if row["aging_deg"] and row["ad_deg"]:
        return "shared"
    if row["aging_deg"]:
        return "aging_only"
    if row["ad_deg"]:
        return "ad_only"
    return "none"


def merge_de(de_aging, de_ad):
    left = de_aging[["gene", "logFC", "deg"]].rename(columns={"logFC": "aging_logFC", "deg": "aging_deg"})
    right = de_ad[["gene", "logFC", "deg"]].rename(columns={"logFC": "ad_logFC", "deg": "ad_deg"})
    merged = left.merge(right, on="gene")
    merged["set"] = merged.apply(assign_set, axis=1)
    same_sign = np.sign(merged["aging_logFC"]) == np.sign(merged["ad_logFC"])
    merged["direction_class"] = np.where(merged["set"] == "shared", np.where(same_sign, "concordant", "discordant"), "")
    return merged


def overlap_stats(merged):
    total = len(merged)
    n_aging = int(merged["aging_deg"].sum())
    n_ad = int(merged["ad_deg"].sum())
    n_shared = int((merged["set"] == "shared").sum())
    union = n_aging + n_ad - n_shared
    rho, rho_p = stats.spearmanr(merged["aging_logFC"], merged["ad_logFC"])
    return {
        "genes_tested": total,
        "aging_degs": n_aging,
        "ad_degs": n_ad,
        "shared": n_shared,
        "aging_only": int((merged["set"] == "aging_only").sum()),
        "ad_only": int((merged["set"] == "ad_only").sum()),
        "concordant": int((merged["direction_class"] == "concordant").sum()),
        "discordant": int((merged["direction_class"] == "discordant").sum()),
        "hypergeom_p": float(stats.hypergeom.sf(n_shared - 1, total, n_aging, n_ad)),
        "jaccard": n_shared / union if union else 0.0,
        "spearman_r": float(rho),
        "spearman_p": float(rho_p),
    }


def ad_only_table(de_ad):
    table = de_ad[de_ad["deg"]][["gene", "logFC"]].rename(columns={"logFC": "ad_logFC"})
    table.insert(1, "aging_logFC", np.nan)
    table["set"] = "ad"
    table["direction_class"] = ""
    return table


def plot_venn(stats_dict, path):
    fig, ax = plt.subplots(figsize=(5, 4))
    venn2(subsets=(stats_dict["aging_only"], stats_dict["ad_only"], stats_dict["shared"]),
          set_labels=("Aging (old vs young)", "AD (AD vs old)"),
          set_colors=(SET_COLORS["aging_only"], SET_COLORS["ad_only"]), ax=ax)
    ax.set_title(f"DEG overlap  (Jaccard {stats_dict['jaccard']:.2f}, p={stats_dict['hypergeom_p']:.1e})")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_scatter(merged, stats_dict, path):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(merged["aging_logFC"], merged["ad_logFC"], s=4, c="#cccccc", linewidths=0, label="not DEG")
    for name, color in SET_COLORS.items():
        sub = merged[merged["set"] == name]
        ax.scatter(sub["aging_logFC"], sub["ad_logFC"], s=8, c=color, linewidths=0, label=f"{name} ({len(sub)})")
    lim = np.nanmax(np.abs(merged[["aging_logFC", "ad_logFC"]].to_numpy())) * 1.05
    ax.plot([-lim, lim], [-lim, lim], color="black", lw=0.6, ls="--")
    ax.axhline(0, color="grey", lw=0.4)
    ax.axvline(0, color="grey", lw=0.4)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("Aging logFC (old vs young)")
    ax.set_ylabel("AD logFC (AD vs old)")
    ax.set_title(f"Spearman r = {stats_dict['spearman_r']:.2f} across {stats_dict['genes_tested']} genes")
    ax.legend(loc="upper left", fontsize=8, markerscale=2)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_compare(de_aging, de_ad, results_dir, figures_dir):
    if de_aging is None:
        table = ad_only_table(de_ad)
        stats_dict = {"genes_tested": len(de_ad), "ad_degs": int(de_ad["deg"].sum()), "aging_available": False}
    else:
        merged = merge_de(de_aging, de_ad)
        stats_dict = overlap_stats(merged)
        stats_dict["aging_available"] = True
        plot_venn(stats_dict, figures_dir / "venn.png")
        plot_scatter(merged, stats_dict, figures_dir / "logfc_scatter.png")
        table = merged[merged["set"] != "none"][["gene", "aging_logFC", "ad_logFC", "set", "direction_class"]]
    table.to_csv(results_dir / "deg_sets.csv", index=False)
    (results_dir / "comparison_stats.json").write_text(json.dumps(stats_dict, indent=2))
    return table, stats_dict


def load_compare(results_dir):
    table = pd.read_csv(results_dir / "deg_sets.csv", keep_default_na=True)
    table["direction_class"] = table["direction_class"].fillna("")
    stats_dict = json.loads((results_dir / "comparison_stats.json").read_text())
    return table, stats_dict
