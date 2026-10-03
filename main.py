import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import config
from src.data import prepare_data, print_groups
from src.de import run_de, mark_degs, plot_volcano
from src.compare import run_compare, load_compare
from src.enrich import run_enrich
from src.network import run_network
from src.report import write_report

ROOT = Path(__file__).parent


def parse_args():
    parser = argparse.ArgumentParser(description="Brain aging vs Alzheimer's disease comparative transcriptomics")
    parser.add_argument("--gse", default=config.GSE)
    parser.add_argument("--region", default=config.REGION)
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--string", action="store_true")
    return parser.parse_args()


def output_dirs(args):
    if args.demo:
        tag = "demo"
    elif args.gse == config.GSE and args.region == config.REGION:
        tag = ""
    else:
        tag = f"{args.gse}_{args.region.replace(' ', '_')}"
    dirs = [ROOT / "results" / tag, ROOT / "figures" / tag, ROOT / "data"]
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)
    return dirs


def demo_data(seed):
    rng = np.random.default_rng(seed)
    genes = [f"GENE{i:04d}" for i in range(2000)]
    groups = ["young"] * 20 + ["old"] * 20 + ["AD"] * 20
    meta = pd.DataFrame({
        "sample": [f"S{i:02d}" for i in range(60)],
        "group": groups,
        "sex": rng.choice(["male", "female"], 60),
        "age": [rng.uniform(20, 39) if g == "young" else rng.uniform(65, 95) for g in groups],
    })
    base = rng.normal(8, 1.5, (2000, 1))
    latent = rng.normal(0, 1, (4, 60))
    values = base + rng.normal(0, 0.3, (2000, 60))
    young, old, ad = (np.array(groups) == g for g in ["young", "old", "AD"])
    plant = {"shared": (0, 120), "discordant": (120, 150), "aging_only": (150, 300), "ad_only": (300, 500)}
    for k, (name, (lo, hi)) in enumerate(plant.items()):
        sign = np.where(np.arange(lo, hi) % 2 == 0, 1.0, -1.0)[:, None]
        values[lo:hi] += 0.25 * sign * latent[k]
        if name in ("shared", "discordant", "aging_only"):
            values[lo:hi][:, old | ad] += 0.8 * sign
        if name == "shared":
            values[lo:hi][:, ad] += 0.8 * sign
        if name == "discordant":
            values[lo:hi][:, ad] -= 1.6 * sign
        if name == "ad_only":
            values[lo:hi][:, ad] += 0.8 * sign
    return pd.DataFrame(values, index=genes, columns=meta["sample"]), meta


def step_data(args, data_dir):
    print(f"[1/6] Loading data: {'synthetic demo' if args.demo else args.gse + ' / ' + args.region}")
    if args.demo:
        expr, meta = demo_data(config.SEED)
        print_groups(meta)
        return expr, meta
    return prepare_data(args.gse, args.region, config, data_dir)


def de_one(expr, meta, name, group_a, group_b, results_dir, figures_dir):
    path = results_dir / f"de_{name}.csv"
    if path.exists():
        return pd.read_csv(path)
    de = mark_degs(run_de(expr, meta, group_a, group_b), config.FDR, config.LOGFC_CUT)
    de.to_csv(path, index=False)
    plot_volcano(de, f"{name}: {group_b} vs {group_a}", figures_dir / f"volcano_{name}.png", config.FDR, config.LOGFC_CUT)
    return de


def step_de(expr, meta, results_dir, figures_dir):
    print("[2/6] Differential expression")
    groups = set(meta["group"])
    de_aging = de_one(expr, meta, "aging", "young", "old", results_dir, figures_dir) if {"young", "old"} <= groups else None
    if de_aging is None:
        print("  no young group: skipping aging comparison")
    de_ad = de_one(expr, meta, "ad", "old", "AD", results_dir, figures_dir)
    print(f"  DEGs: aging {int(de_aging['deg'].sum()) if de_aging is not None else 'n/a'}, AD {int(de_ad['deg'].sum())}")
    return de_aging, de_ad


def step_compare(de_aging, de_ad, results_dir, figures_dir):
    print("[3/6] Comparing aging vs AD")
    if (results_dir / "deg_sets.csv").exists() and (results_dir / "comparison_stats.json").exists():
        return load_compare(results_dir)
    return run_compare(de_aging, de_ad, results_dir, figures_dir)


def step_enrich(table, set_names, demo, results_dir, figures_dir):
    print("[4/6] Pathway enrichment" + (" (skipped in demo mode)" if demo else ""))
    results = {}
    for name in set_names:
        path = results_dir / f"enrichment_{name}.csv"
        if path.exists() and not demo:
            results[name] = pd.read_csv(path)
        else:
            results[name] = run_enrich(table, name, config.GENE_SETS, results_dir, figures_dir, skip=demo)
    return results


def step_network(expr, meta, table, set_names, use_string, results_dir, figures_dir):
    print("[5/6] Co-expression hub genes")
    path = results_dir / "hubs.csv"
    if path.exists() and (not use_string or "string_degree" in pd.read_csv(path, nrows=0).columns):
        return pd.read_csv(path)
    hubs = run_network(expr, meta, table, set_names, config, figures_dir, use_string)
    hubs.to_csv(path, index=False)
    return hubs


def main():
    args = parse_args()
    results_dir, figures_dir, data_dir = output_dirs(args)
    expr, meta = step_data(args, data_dir)
    de_aging, de_ad = step_de(expr, meta, results_dir, figures_dir)
    table, stats = step_compare(de_aging, de_ad, results_dir, figures_dir)
    set_names = ["shared", "aging_only", "ad_only"] if stats["aging_available"] else ["ad"]
    enrichments = step_enrich(table, set_names, args.demo, results_dir, figures_dir)
    hubs = step_network(expr, meta, table, set_names, args.string, results_dir, figures_dir)
    print("[6/6] Writing report")
    write_report(results_dir / "report.md", "synthetic demo" if args.demo else args.gse, args.region, meta, stats, enrichments, hubs)
    print(f"Done. Results in {results_dir.relative_to(ROOT)}, figures in {figures_dir.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
