#!/usr/bin/env python3
"""Step 2. Summarise raw and corrected branch variation; derive simulation parameters.

Outputs (in --out):
  species_terminal_stats.csv   per-species raw (n, mean, SD, median, IQR, CV) + corrected residual stats (sigma_s)
  pooled_sigma.csv             sigma_external, sigma_internal, sigma_gene (+ MAD-based, + gene-cluster bootstrap CI)
  sigma_by_length_bin.csv      heteroscedasticity check (does sigma depend on reference branch length?)
  sigma_by_clade_size.csv      same for internal branches vs clade size
  sim_params.json              input for 03_simulate_trees.py
  diagnostics.png
"""
import argparse, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

K = 1.4826  # MAD -> SD under normality


def sd(x):
    x = np.asarray(x, float); x = x[~np.isnan(x)]
    return float(np.std(x, ddof=1)) if len(x) > 1 else np.nan


def mad_sd(x):
    x = np.asarray(x, float); x = x[~np.isnan(x)]
    return float(K * np.median(np.abs(x - np.median(x)))) if len(x) > 1 else np.nan


def iqr(x):
    x = np.asarray(x, float); x = x[~np.isnan(x)]
    return float(np.percentile(x, 75) - np.percentile(x, 25)) if len(x) else np.nan


def species_table(ter, min_n):
    rows = []
    for sp, d in ter.groupby("species"):
        L = d["L_gene"].values
        E = d["E"].dropna().values
        mean = L.mean()
        s = sd(L)
        rows.append(dict(species=sp, n_gene_trees=len(L), mean_L=mean, sd_L=s, median_L=np.median(L),
                         iqr_L=iqr(L), cv_L=s / mean if mean > 0 else np.nan,
                         n_E=len(E), mean_E=E.mean() if len(E) else np.nan,
                         sigma_s=sd(E), sigma_s_mad=mad_sd(E), median_E=np.median(E) if len(E) else np.nan,
                         iqr_E=iqr(E), low_n=int(len(E) < min_n)))
    return pd.DataFrame(rows).sort_values("sigma_s", ascending=False)


def binned(df, col, bins, label):
    d = df.dropna(subset=["E"]).copy()
    if col == "L_ref":
        d["bin"] = pd.qcut(np.log10(d[col]), bins, duplicates="drop")
    else:
        d["bin"] = pd.qcut(d[col], bins, duplicates="drop")
    out = d.groupby("bin", observed=True).agg(n=("E", "size"), sd_E=("E", sd), mad_sd_E=("E", mad_sd),
                                              median_ref=(col, "median")).reset_index()
    out["bin"] = out["bin"].astype(str); out.insert(0, "branch_class", label)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-n-species", type=int, default=10, help="species with fewer observations fall back to pooled sigma")
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    import os; os.makedirs(a.out, exist_ok=True)

    ter = pd.read_csv(f"{a.in_dir}/terminal_branches.csv")
    inn = pd.read_csv(f"{a.in_dir}/internal_branches.csv")
    gs = pd.read_csv(f"{a.in_dir}/gene_summary.csv")

    sp = species_table(ter, a.min_n_species)
    sp.to_csv(f"{a.out}/species_terminal_stats.csv", index=False)

    Et, Ei = ter["E"].dropna().values, inn["E"].dropna().values
    Mg = gs["M_g"].dropna().values
    # variance without species-specific offsets (how much of the pooled SD is lineage-level bias?)
    Et_dm = ter.dropna(subset=["E"]).assign(Ed=lambda d: d["E"] - d.groupby("species")["E"].transform("mean"))["Ed"].values

    # gene-cluster bootstrap (residuals within a gene are not independent)
    rng = np.random.default_rng(a.seed)
    tg = {g: d["E"].dropna().values for g, d in ter.groupby("gene")}
    ig = {g: d["E"].dropna().values for g, d in inn.groupby("gene")}
    genes = list(gs["gene"])
    bt, bi, bg = [], [], []
    for _ in range(a.n_boot):
        pick = rng.choice(len(genes), len(genes), replace=True)
        bt.append(sd(np.concatenate([tg.get(genes[i], np.array([])) for i in pick])))
        bi.append(sd(np.concatenate([ig.get(genes[i], np.array([])) for i in pick])))
        bg.append(sd(gs["M_g"].values[pick]))
    ci = lambda b: (float(np.nanpercentile(b, 2.5)), float(np.nanpercentile(b, 97.5)))

    pooled = pd.DataFrame([
        dict(parameter="sigma_external", n=len(Et), sd=sd(Et), mad_sd=mad_sd(Et), ci_lo=ci(bt)[0], ci_hi=ci(bt)[1],
             mean=float(Et.mean())),
        dict(parameter="sigma_external_species_demeaned", n=len(Et_dm), sd=sd(Et_dm), mad_sd=mad_sd(Et_dm),
             ci_lo=np.nan, ci_hi=np.nan, mean=0.0),
        dict(parameter="sigma_internal", n=len(Ei), sd=sd(Ei), mad_sd=mad_sd(Ei), ci_lo=ci(bi)[0], ci_hi=ci(bi)[1],
             mean=float(Ei.mean()) if len(Ei) else np.nan),
        dict(parameter="sigma_gene (SD of M_g)", n=len(Mg), sd=sd(Mg), mad_sd=mad_sd(Mg), ci_lo=ci(bg)[0], ci_hi=ci(bg)[1],
             mean=float(Mg.mean())),
    ])
    pooled.to_csv(f"{a.out}/pooled_sigma.csv", index=False)

    pd.concat([binned(ter, "L_ref", 5, "terminal"), binned(inn, "L_ref", 5, "internal")]).to_csv(
        f"{a.out}/sigma_by_length_bin.csv", index=False)
    binned(inn, "clade_size", 5, "internal").to_csv(f"{a.out}/sigma_by_clade_size.csv", index=False)

    params = dict(
        sigma_external=sd(Et), sigma_internal=sd(Ei), sigma_gene=sd(Mg), mu_M=float(Mg.mean()),
        sigma_external_robust=mad_sd(Et), sigma_internal_robust=mad_sd(Ei), sigma_gene_robust=mad_sd(Mg),
        species_min_n=a.min_n_species,
        species_sigma={r.species: r.sigma_s for r in sp.itertuples() if not r.low_n and not np.isnan(r.sigma_s)},
        species_sigma_robust={r.species: r.sigma_s_mad for r in sp.itertuples() if not r.low_n and not np.isnan(r.sigma_s_mad)},
    )
    json.dump(params, open(f"{a.out}/sim_params.json", "w"), indent=2)

    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    ax[0, 0].hist(Et, bins=60, color="#3b6ea5"); ax[0, 0].set_title(f"Terminal E (SD={sd(Et):.3f}, MADsd={mad_sd(Et):.3f})")
    ax[0, 1].hist(Ei, bins=60, color="#c0572b"); ax[0, 1].set_title(f"Internal E (SD={sd(Ei):.3f}, MADsd={mad_sd(Ei):.3f})")
    ax[0, 2].hist(Mg, bins=40, color="#666"); ax[0, 2].set_title(f"Gene-wide M_g (SD={sd(Mg):.3f})")
    d = ter.dropna(subset=["E"]); ax[1, 0].scatter(np.log10(d["L_ref"]), d["E"], s=3, alpha=.3, color="#3b6ea5")
    ax[1, 0].set_xlabel("log10 reference terminal length"); ax[1, 0].set_ylabel("E")
    d = inn.dropna(subset=["E"]); ax[1, 1].scatter(np.log10(d["L_ref"]), d["E"], s=3, alpha=.3, color="#c0572b")
    ax[1, 1].set_xlabel("log10 reference internal length"); ax[1, 1].set_ylabel("E")
    s2 = sp.dropna(subset=["sigma_s"]); ax[1, 2].bar(range(len(s2)), s2["sigma_s"], color=np.where(s2["low_n"] == 1, "#bbb", "#3b6ea5"))
    ax[1, 2].axhline(sd(Et), color="k", ls="--"); ax[1, 2].set_title("sigma_s per species (grey: n too small)")
    plt.tight_layout(); plt.savefig(f"{a.out}/diagnostics.png", dpi=150)
    print(pooled.to_string(index=False))


if __name__ == "__main__":
    main()
