import pandas as pd

LIMITATIONS = [
    "Single brain region per run.",
    "Single dataset; results need replication in an independent cohort.",
    "\"Early AD\" depends on how the dataset labels its AD cases; many are late-stage post-mortem diagnoses.",
    "Post-mortem tissue: agonal state, RNA quality and post-mortem interval are not modeled.",
    "Cell-type composition changes (neuron loss, gliosis) are not modeled and can drive expression differences.",
    "Co-expression hubs are correlational, not causal.",
    "Batch is not modeled; in GSE48350 all AD samples were submitted separately from the controls, so AD vs old is confounded with batch.",
]


def interpret(s):
    if not s.get("aging_available"):
        return "Aging comparison not available (no young controls in this dataset); only AD vs old controls is reported."
    concordant_frac = s["concordant"] / s["shared"] if s["shared"] else 0.0
    discordant_frac = s["discordant"] / s["shared"] if s["shared"] else 0.0
    if (s["jaccard"] >= 0.3 or s["spearman_r"] >= 0.5) and concordant_frac >= 0.8:
        return "AD looks largely like accelerated aging in this region."
    if s["jaccard"] < 0.15 and s["ad_only"] > 2 * s["shared"] and (s["shared"] == 0 or discordant_frac >= 0.2):
        return "AD shows a distinct molecular signature beyond aging."
    return "Mixed: partial overlap with substantial AD-specific changes."


def group_lines(meta):
    counts = meta["group"].value_counts()
    return [f"- {g}: {int(counts.get(g, 0))}" for g in ["young", "old", "AD"]]


def stats_lines(s):
    if not s.get("aging_available"):
        return [f"- AD DEGs: {s['ad_degs']} of {s['genes_tested']} genes tested"]
    return [
        f"- Genes tested in both: {s['genes_tested']}",
        f"- Aging DEGs: {s['aging_degs']}",
        f"- AD DEGs: {s['ad_degs']}",
        f"- Shared: {s['shared']} (concordant {s['concordant']}, discordant {s['discordant']})",
        f"- Aging-only: {s['aging_only']}",
        f"- AD-only: {s['ad_only']}",
        "",
        f"- Overlap hypergeometric p: {s['hypergeom_p']:.2e}",
        f"- Jaccard index: {s['jaccard']:.3f}",
        f"- Spearman r of logFC (all common genes): {s['spearman_r']:.3f}",
    ]


def pathway_lines(enrichments):
    lines = []
    for set_name, result in enrichments.items():
        lines.append(f"**{set_name}**")
        if result.empty:
            lines.append("- no enrichment results")
        for _, row in result.head(5).iterrows():
            lines.append(f"- {row['term']} ({row['library']}, {row['direction']}, adj p {row['adj_p']:.1e})")
        lines.append("")
    return lines


def hub_lines(hubs, set_names):
    lines = []
    for set_name in set_names:
        sub = hubs[hubs["set"] == set_name]
        top = ", ".join(f"{g} ({d})" for g, d in zip(sub["gene"].head(5), sub["degree"].head(5)))
        lines.append(f"- **{set_name}**: {top or 'none'}")
    return lines


def write_report(path, gse, region, meta, s, enrichments, hubs):
    lines = [
        "# Brain aging vs Alzheimer's disease: summary",
        "",
        f"Dataset: {gse}  |  Region: {region}",
        "",
        "## Group sizes",
        *group_lines(meta),
        "",
        "## Differential expression and overlap",
        *stats_lines(s),
        "",
        "## Top pathways",
        *pathway_lines(enrichments),
        "## Top hubs (co-expression degree)",
        *hub_lines(hubs, list(enrichments)),
        "",
        "## Reading the result",
        interpret(s),
        "",
        "## Limitations",
        *[f"- {x}" for x in LIMITATIONS],
        "",
    ]
    path.write_text("\n".join(lines))
