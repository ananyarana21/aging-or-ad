import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests


def has_sex(meta):
    return meta["sex"].notna().all() and meta["sex"].nunique() == 2


def ols_test(values, meta, group_b):
    group = (meta["group"] == group_b).astype(float).to_numpy()
    sex = (meta["sex"] == meta["sex"].iloc[0]).astype(float).to_numpy()
    x = np.column_stack([np.ones(len(group)), group, sex])
    xtx_inv = np.linalg.inv(x.T @ x)
    with np.errstate(all="ignore"):
        beta = xtx_inv @ x.T @ values.T
        resid = values.T - x @ beta
    df = x.shape[0] - x.shape[1]
    sigma2 = (resid ** 2).sum(axis=0) / df
    se = np.sqrt(sigma2 * xtx_inv[1, 1])
    t = beta[1] / se
    p = 2 * stats.t.sf(np.abs(t), df)
    return t, p


def welch_test(values_a, values_b):
    t, p = stats.ttest_ind(values_b, values_a, axis=1, equal_var=False)
    return t, p


def run_de(expr, meta, group_a, group_b):
    meta = meta[meta["group"].isin([group_a, group_b])]
    values = expr[meta["sample"]].to_numpy()
    in_a = (meta["group"] == group_a).to_numpy()
    mean_a = values[:, in_a].mean(axis=1)
    mean_b = values[:, ~in_a].mean(axis=1)
    if has_sex(meta):
        t, p = ols_test(values, meta, group_b)
    else:
        t, p = welch_test(values[:, in_a], values[:, ~in_a])
    p = np.nan_to_num(p, nan=1.0)
    padj = multipletests(p, method="fdr_bh")[1]
    return pd.DataFrame({
        "gene": expr.index,
        "logFC": mean_b - mean_a,
        "mean_a": mean_a,
        "mean_b": mean_b,
        "t": t,
        "pval": p,
        "padj": padj,
    })


def mark_degs(de, fdr, logfc_cut):
    de["deg"] = (de["padj"] < fdr) & (de["logFC"].abs() > logfc_cut)
    return de


def plot_volcano(de, title, path, fdr, logfc_cut):
    y = -np.log10(de["padj"].clip(lower=1e-300))
    colors = np.where(de["deg"], np.where(de["logFC"] > 0, "#c0392b", "#2c7bb6"), "#bbbbbb")
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(de["logFC"], y, c=colors, s=6, linewidths=0)
    ax.axhline(-np.log10(fdr), color="black", lw=0.6, ls="--")
    ax.axvline(logfc_cut, color="black", lw=0.6, ls="--")
    ax.axvline(-logfc_cut, color="black", lw=0.6, ls="--")
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p")
    ax.set_title(f"{title}  ({int(de['deg'].sum())} DEGs)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
