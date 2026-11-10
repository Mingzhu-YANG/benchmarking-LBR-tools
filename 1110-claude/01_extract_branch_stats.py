#!/usr/bin/env python3
"""Step 1. Extract terminal and split-matched internal branch lengths from gene trees,
compare them with the reference tree pruned to each gene's taxon set, estimate the
gene-wide rate factor M_g and compute corrected residuals E = log(L/Lref) - M_g.

Usage
  python 01_extract_branch_stats.py --ref reference.tre --genes gene_trees/ --out results/
  python 01_extract_branch_stats.py --ref reference.tre --genes all_genes.tre --out results/   # one tree per line
"""
import argparse, glob, os, sys
import numpy as np
import pandas as pd
from bl_utils import read_newick, raw_edges, restrict, mask_to_labels, popcount


def load_gene_trees(path, pattern, strip_re):
    items = []
    if os.path.isdir(path):
        for f in sorted(glob.glob(os.path.join(path, pattern))):
            gid = os.path.basename(f).split(".")[0]
            items.append((gid, read_newick(open(f).read(), strip_re)))
    else:
        lines = [l.strip() for l in open(path) if l.strip()]
        for i, l in enumerate(lines, 1):
            items.append((f"g{i:04d}", read_newick(l, strip_re)))
    return items


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--genes", required=True, help="directory of gene tree files, or one file with one newick per line")
    ap.add_argument("--pattern", default="*", help="glob inside --genes directory (e.g. '*.treefile')")
    ap.add_argument("--out", required=True)
    ap.add_argument("--strip-regex", default=None, help="regex removed from tip labels (e.g. '_[0-9]+$')")
    ap.add_argument("--min-length", type=float, default=1e-6,
                    help="branches shorter than this (either tree) are censored: no log ratio, counted separately")
    ap.add_argument("--m-method", choices=["median", "logtl"], default="median",
                    help="gene-wide factor: median of log ratios over all matched branches (default) or log of total-length ratio")
    ap.add_argument("--min-branches-for-M", type=int, default=10)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    ref = read_newick(open(a.ref).read(), a.strip_regex)
    labels = sorted(lf.taxon.label for lf in ref.leaf_node_iter())
    index = {l: i for i, l in enumerate(labels)}
    ref_edges, _ = raw_edges(ref, index)
    genes = load_gene_trees(a.genes, a.pattern, a.strip_regex)
    print(f"reference: {len(labels)} taxa; gene trees: {len(genes)}", file=sys.stderr)

    ter_rows, int_rows, g_rows, taxa_rows = [], [], [], []
    for gid, gt in genes:
        tips = [lf.taxon.label for lf in gt.leaf_node_iter()]
        missing = [t for t in tips if t not in index]
        if missing:
            print(f"[warn] {gid}: {len(missing)} tips not in reference ignored (e.g. {missing[:3]})", file=sys.stderr)
        g_edges, S = raw_edges(gt, index)
        n = popcount(S)
        if n < 4:
            print(f"[skip] {gid}: <4 usable taxa", file=sys.stderr); continue
        g_term, g_int = restrict(g_edges, S, labels)     # gene tree (unrooted, split-keyed)
        r_term, r_int = restrict(ref_edges, S, labels)   # reference pruned to S

        matched = sorted(set(g_int) & set(r_int))
        ref_only = len(set(r_int) - set(g_int))          # topology mismatch, NOT zero length
        gene_only = len(set(g_int) - set(r_int))

        rows_t, rows_i = [], []
        for sp, L in g_term.items():
            Lr = r_term.get(sp)
            ok = Lr is not None and L > a.min_length and Lr > a.min_length
            rows_t.append(dict(gene=gid, species=sp, L_gene=L, L_ref=Lr,
                               log_ratio=np.log(L / Lr) if ok else np.nan, censored=int(not ok)))
        for m in matched:
            L, Lr = g_int[m], r_int[m]
            ok = L > a.min_length and Lr > a.min_length
            small = mask_to_labels(m, labels)
            rows_i.append(dict(gene=gid, split_id=hex(m), clade_size=min(len(small), n - len(small)),
                               L_gene=L, L_ref=Lr,
                               log_ratio=np.log(L / Lr) if ok else np.nan, censored=int(not ok)))

        lr = np.array([r["log_ratio"] for r in rows_t + rows_i], float)
        lr = lr[~np.isnan(lr)]
        Ls = sum(r["L_gene"] for r in rows_t + rows_i if r["censored"] == 0)
        Lrs = sum(r["L_ref"] for r in rows_t + rows_i if r["censored"] == 0)
        M_med = float(np.median(lr)) if len(lr) >= a.min_branches_for_M else np.nan
        M_tl = float(np.log(Ls / Lrs)) if len(lr) >= a.min_branches_for_M and Lrs > 0 else np.nan
        M = M_med if a.m_method == "median" else M_tl
        for r in rows_t + rows_i:
            r["M_g"] = M
            r["E"] = r["log_ratio"] - M
        ter_rows += rows_t; int_rows += rows_i
        n_ref_int = len(r_int)
        g_rows.append(dict(gene=gid, n_taxa=n, n_ref_internal=n_ref_int, n_gene_internal=len(g_int),
                           n_matched_internal=len(matched), frac_ref_splits_recovered=len(matched) / max(n_ref_int, 1),
                           n_ref_only_splits=ref_only, n_gene_only_splits=gene_only,
                           n_censored_terminal=sum(r["censored"] for r in rows_t),
                           n_censored_internal=sum(r["censored"] for r in rows_i),
                           n_branches_for_M=len(lr), M_median=M_med, M_logTL=M_tl, M_g=M))
        taxa_rows.append((gid, ",".join(mask_to_labels(S, labels))))

    pd.DataFrame(ter_rows).to_csv(f"{a.out}/terminal_branches.csv", index=False)
    pd.DataFrame(int_rows).to_csv(f"{a.out}/internal_branches.csv", index=False)
    pd.DataFrame(g_rows).to_csv(f"{a.out}/gene_summary.csv", index=False)
    with open(f"{a.out}/gene_taxa.tsv", "w") as f:
        for gid, t in taxa_rows:
            f.write(f"{gid}\t{t}\n")
    print(f"wrote tables to {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
