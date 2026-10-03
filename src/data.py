import re
import numpy as np
import pandas as pd
import GEOparse

KEY_PATTERNS = {
    "region": re.compile(r"\bregion\b", re.I),
    "disease": re.compile(r"disease|diagnosis|status|condition", re.I),
    "age": re.compile(r"\bage\b", re.I),
    "sex": re.compile(r"\bsex\b|gender", re.I),
}
AD_PATTERN = re.compile(r"(?i:alzheimer)|(?:^|[\s_,;:(-])AD(?:[\s_,;:)-]|$)")
CONTROL_PATTERN = re.compile(r"control|normal|healthy", re.I)
AGE_PATTERN = re.compile(r"(\d+(?:\.\d+)?)")
SEX_PATTERN = re.compile(r"^\s*(f|female|m|male)\b", re.I)


def load_gse(gse, data_dir):
    return GEOparse.get_GEO(geo=gse, destdir=str(data_dir), silent=True)


def parse_characteristics(entries):
    fields = {}
    for entry in entries:
        key, _, value = entry.partition(":")
        for field, pattern in KEY_PATTERNS.items():
            if field not in fields and pattern.search(key):
                fields[field] = value.strip()
    return fields


def sample_record(name, gsm):
    fields = parse_characteristics(gsm.metadata.get("characteristics_ch1", []))
    title = gsm.metadata.get("title", [""])[0]
    return {
        "sample": name,
        "title": title,
        "region_raw": fields.get("region", gsm.metadata.get("source_name_ch1", [""])[0]),
        "disease_raw": fields.get("disease", ""),
        "age_raw": fields.get("age", ""),
        "sex_raw": fields.get("sex", ""),
    }


def disease_state(row):
    text = row["disease_raw"] or row["title"]
    if AD_PATTERN.search(text):
        return "AD"
    if row["disease_raw"] and not CONTROL_PATTERN.search(row["disease_raw"]):
        return "other"
    return "control"


def parse_age(raw):
    match = AGE_PATTERN.search(raw)
    return float(match.group(1)) if match else np.nan


def parse_sex(raw):
    match = SEX_PATTERN.search(raw)
    return ("female" if match.group(1).lower().startswith("f") else "male") if match else np.nan


def build_meta(geo):
    meta = pd.DataFrame([sample_record(n, g) for n, g in geo.gsms.items()])
    meta["region"] = meta["region_raw"].str.strip().str.lower()
    meta["disease"] = meta.apply(disease_state, axis=1)
    meta["age"] = meta["age_raw"].map(parse_age)
    meta["sex"] = meta["sex_raw"].map(parse_sex)
    return meta


def print_raw_values(meta):
    for col in ["region_raw", "disease_raw", "sex_raw"]:
        print(f"  {col}: {sorted(meta[col].unique().tolist())}")
    ages = meta["age"].dropna()
    print(f"  age_raw examples: {meta['age_raw'].unique()[:8].tolist()} -> parsed {len(ages)}/{len(meta)}, range {ages.min():.0f}-{ages.max():.0f}")
    print(f"  disease (title fallback used where no field): {meta['disease'].value_counts().to_dict()}")


def assign_group(row, young_max, old_min):
    if row["disease"] == "AD":
        return "AD"
    if row["disease"] != "control" or np.isnan(row["age"]):
        return None
    if row["age"] < young_max:
        return "young"
    if row["age"] >= old_min:
        return "old"
    return None


def probe_symbols(geo):
    gpl = next(iter(geo.gpls.values())).table
    col = next(c for c in gpl.columns if re.fullmatch(r"gene[ _]?symbol", c, re.I))
    symbols = gpl.set_index("ID")[col].dropna().astype(str).str.strip()
    return symbols[(symbols != "") & ~symbols.str.contains("///")]


def build_expr(geo, samples):
    expr = geo.pivot_samples("VALUE")[samples].apply(pd.to_numeric, errors="coerce").dropna()
    if expr.to_numpy().max() > 100:
        expr = np.log2(expr.clip(lower=0) + 1)
    symbols = probe_symbols(geo)
    expr = expr.loc[expr.index.intersection(symbols.index)]
    expr["gene"] = symbols.loc[expr.index]
    expr["mean"] = expr[samples].mean(axis=1)
    expr = expr.sort_values("mean", ascending=False).drop_duplicates("gene")
    return expr.set_index("gene")[samples]


def print_groups(meta):
    counts = meta["group"].value_counts()
    print(f"  group sizes: {counts.to_dict()}")
    for group in ["young", "old", "AD"]:
        if counts.get(group, 0) < 5:
            print(f"  WARNING: group '{group}' has only {counts.get(group, 0)} samples")


def prepare_data(gse, region, cfg, data_dir):
    expr_path = data_dir / f"expr_{gse}_{region}.csv"
    meta_path = data_dir / f"meta_{gse}_{region}.csv"
    if expr_path.exists() and meta_path.exists():
        expr, meta = pd.read_csv(expr_path, index_col=0), pd.read_csv(meta_path)
    else:
        geo = load_gse(gse, data_dir)
        meta = build_meta(geo)
        print_raw_values(meta)
        meta = meta[meta["region"].str.contains(region.lower(), regex=False)].copy()
        meta["group"] = meta.apply(assign_group, axis=1, args=(cfg.YOUNG_MAX_AGE, cfg.OLD_MIN_AGE))
        meta = meta.dropna(subset=["group"]).reset_index(drop=True)
        expr = build_expr(geo, meta["sample"].tolist())
        expr.to_csv(expr_path)
        meta.to_csv(meta_path, index=False)
    print_groups(meta)
    return expr, meta
