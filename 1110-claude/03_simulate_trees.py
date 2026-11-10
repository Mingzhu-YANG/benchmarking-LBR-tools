#!/usr/bin/env python3
"""Step 3. Simulate gene-tree branch lengths by perturbing the (pruned) reference tree.

  L_sim = L_ref * exp( M_g + eps ),   M_g ~ N(mu_M, sigma_gene),
  eps ~ N(0, sigma_external) on terminal branches (or sigma_s per species with --mode species),
  eps ~ N(0, sigma_internal) on internal branches.

--taxa-sets gene_taxa.tsv (from step 1) makes each replicate reuse a real gene's taxon set
(sampled with replacement), so replicates have 30-60 taxa like the empirical data.
--scale multiplies all sigmas (sensitivity analysis).  --robust uses MAD-based sigmas.
--mean-preserving subtracts sigma^2/2 so that E[L_sim] = L_ref * exp(M_g).
"""
import argparse, json, os
import numpy as np
import pandas as pd
import dendropy
from bl_utils import read_newick


def unroot_sum(tree):
    """Make the root trifurcating; the two root edges are merged into one (lengths summed)."""
    if len(tree.seed_node.child_nodes()) == 2:
        tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    tree.is_rooted = False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--out", required=True)
    ap.add_argument("--taxa-sets", default=None)
    ap.add_argument("--mode", choices=["pooled", "species"], default="pooled")
    ap.add_argument("--robust", action="store_true")
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--mean-preserving", action="store_true")
    ap.add_argument("--no-gene-rate", action="store_true", help="set M_g = 0")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    P = json.load(open(a.params))
    r = "_robust" if a.robust else ""
    s_ext, s_int, s_gene = (P[f"sigma_external{r}"] * a.scale, P[f"sigma_internal{r}"] * a.scale, P[f"sigma_gene{r}"])
    s_sp = {k: v * a.scale for k, v in P[f"species_sigma{r}"].items()}
    mu_M = P["mu_M"]

    ref = read_newick(open(a.ref).read())
    all_taxa = [lf.taxon.label for lf in ref.leaf_node_iter()]
    sets = None
    if a.taxa_sets:
        sets = [l.rstrip("\n").split("\t")[1].split(",") for l in open(a.taxa_sets) if l.strip()]

    truth = []
    for i in range(1, a.n + 1):
        S = [t for t in (sets[rng.integers(len(sets))] if sets else all_taxa) if t in set(all_taxa)]
        t = ref.clone(depth=2)
        t.retain_taxa_with_labels(S)          # prunes and sums lengths through suppressed degree-2 nodes
        unroot_sum(t)
        M = 0.0 if a.no_gene_rate else rng.normal(mu_M, s_gene)
        for nd in t.postorder_node_iter():
            if nd is t.seed_node or nd.edge_length is None:
                continue
            if nd.is_leaf():
                s = s_sp.get(nd.taxon.label, s_ext) if a.mode == "species" else s_ext
            else:
                s = s_int
            eps = rng.normal(0.0, s) - (s * s / 2 if a.mean_preserving else 0.0)
            nd.edge_length = nd.edge_length * float(np.exp(M + eps))
        for nd in t.preorder_internal_node_iter():
            nd.label = None
        with open(f"{a.out}/sim_{i:03d}.tre", "w") as f:
            f.write(t.as_string("newick", suppress_rooting=True))
        truth.append(dict(rep=i, n_taxa=len(S), M_true=M))
    pd.DataFrame(truth).to_csv(f"{a.out}/truth.csv", index=False)
    print(f"wrote {a.n} trees to {a.out}")


if __name__ == "__main__":
    main()
