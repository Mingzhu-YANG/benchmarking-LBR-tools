# Calibrating Branch-Length Variation for Gene-Tree Simulation from Empirical Gene Trees

## 1. Objective

Simulate 100 gene trees by perturbing the branch lengths of a reference species tree (`reference.tre`), with perturbation magnitudes estimated from 301 empirical gene trees (30–60 taxa each). The empirical data are used to answer three questions:

1. How variable are **terminal (external) branches**, pooled and per species? (`σ_external`, `σ_s`)
2. How variable are **internal branches** that are present in both the gene tree and the reference? (`σ_internal`)
3. How variable is the **gene-wide rate** itself? (`σ_gene`, the SD of `M_g`). This is the component removed by the correction and must be put back during simulation.

Topology is out of scope here (discordance is generated elsewhere, e.g. by the coalescent). Topological mismatch is only *reported* and never coded as a zero-length branch.

## 2. Inputs and assumptions

| Item | Requirement |
|---|---|
| `reference.tre` | Newick with branch lengths, tip labels covering all gene-tree taxa. Rooted or unrooted. |
| Gene trees | 301 trees, 30–60 taxa, branch lengths in the **same units** as the reference (substitutions/site). Directory of files or one newick per line. |
| Units | If the reference is in other units (e.g. coalescent units), `M_g` absorbs only a *global* rescaling, so lineage-specific differences will leak into the residuals. Rescale or re-estimate reference branch lengths on the same data first. |
| Tip labels | Must match between trees (use `--strip-regex` for allele/copy suffixes). |

## 3. Definitions

Let `L_{g,b}` be the length of branch `b` in gene tree `g`, and `L_ref,g,b` the length of the corresponding branch in the reference **pruned to the taxa of gene `g`**.

**Matching.** Everything is unrooted and matched by taxon bipartition (split), never by node number.

- Pruning is implemented as projecting every reference edge onto the gene's taxon set. Edges that induce the same bipartition (taxa pruned in between, or the two root edges of a rooted tree) are **summed**, which equals pruning + suppressing degree-2 nodes + unrooting.
- Terminal branch: split with one taxon on one side. Internal branch: all other splits.
- Only splits present in **both** trees enter the internal-branch statistics. A reference split absent from the gene tree is a topological mismatch and is counted (`n_ref_only_splits`) but never assigned length 0.

**Log ratio:** `R_{g,b} = log( L_{g,b} / L_ref,g,b )`

**Gene-wide factor (default):** `M_g = median_b R_{g,b}` over all matched terminal + internal branches. Sensitivity alternative: `M_g = log( Σ L_{g,b} / Σ L_ref,g,b )`.

**Corrected residual:**

$$E_{g,b} = \log\left(\frac{L_{g,b}}{L_{\mathrm{ref},g,b}}\right) - M_g$$

**Per-species raw statistics** (terminal branches, uncorrected): `n_s`, mean `L̄_s`, SD, median, IQR, CV

$$SD_s=\sqrt{\frac{\sum_{g=1}^{n_s}(L_{g,s}-\bar L_s)^2}{n_s-1}},\qquad CV_s=\frac{SD_s}{\bar L_s}$$

**Per-species corrected SD:** `σ_s = SD({E_{g,s}}_g)` (species with `n_s < 10` fall back to the pooled value).

**Pooled parameters:**

| Parameter | Definition | Use |
|---|---|---|
| `σ_external` | SD of all terminal `E_{g,s}` | external-branch noise |
| `σ_internal` | SD of all matched internal `E_{g,b}` | internal-branch noise |
| `σ_s` | SD of `E_{g,s}` per species | species-specific terminal noise / which taxa are unstable |
| `σ_gene` (+ `μ_M`) | SD (and mean) of `M_g` across genes | gene-wide rate variation |
| MAD-based versions | `1.4826 · median(|E − median(E)|)` | robustness check |

## 4. Workflow

```
reference.tre + 301 gene trees
   │  01_extract_branch_stats.py   prune/project ref, split-match, M_g, E
   ▼
terminal_branches.csv  internal_branches.csv  gene_summary.csv  gene_taxa.tsv
   │  02_summarise_sigma.py        raw stats, sigma_s, pooled sigma, MAD, bootstrap CI, diagnostics
   ▼
species_terminal_stats.csv  pooled_sigma.csv  sigma_by_*.csv  sim_params.json  diagnostics.png
   │  (inspect, decide parameters — Section 5)
   │  03_simulate_trees.py         L_sim = L_ref · exp(M_g + ε), 100 trees
   ▼
results/sim_*/                                   
   │  01 + 02 again on simulated trees           self-consistency check
   ▼
estimated σ ≈ input σ ?
```

Run everything: `bash run_all.sh reference.tre gene_trees/ "*.treefile"`.

### Step 1: extraction (`01_extract_branch_stats.py`)
- Parse reference and gene trees; compute edge-to-taxon-set masks once.
- For each gene: project reference onto the gene's taxa, extract terminal and internal edges from both trees, intersect internal splits.
- Branches with length `< --min-length` (default 1e-6) in either tree are **censored** (no log ratio, counted per gene). Zero-length internal branches are typically collapsed/unresolved nodes, not true zeros.
- Estimate `M_g` (needs ≥ 10 usable branches), compute `E` for every branch.
- Per-gene report: taxa, reference splits recovered (`frac_ref_splits_recovered`), mismatches, censored counts, `M_g`.

### Step 2: summaries (`02_summarise_sigma.py`)
- Per-species table: raw (`n`, mean, SD, median, IQR, CV) and corrected (`mean_E`, `σ_s`, MAD-`σ_s`, IQR of `E`).
- Pooled `σ_external`, `σ_internal`, `σ_gene` with **gene-cluster bootstrap** 95% CIs (residuals within a gene are not independent).
- `σ_external` after removing species means (`species_demeaned`): the gap between this and `σ_external` is the share of variance that is a stable lineage effect rather than gene-to-gene noise.
- Heteroscedasticity: `σ` by reference-length quintile and by clade size.

### Step 3: simulation (`03_simulate_trees.py`)
For each of 100 replicates: sample a real gene's taxon set (keeps 30–60 taxa), prune the reference, draw `M_g ~ N(μ_M, σ_gene)`, then multiply each branch by `exp(M_g + ε)` with `ε ~ N(0, σ_external or σ_s)` for terminals and `N(0, σ_internal)` for internals. The two root edges are merged first so one unrooted edge gets one draw. Options: `--mode species`, `--robust`, `--scale k`, `--mean-preserving`, `--no-gene-rate`.

## 5. How to decide the final σ values

1. **Report all three blocks** (`σ_external`, `σ_internal`, `σ_s`) with CIs, plus `σ_gene`.
2. **Compare SD with MAD-based σ.** If SD ≫ MAD-σ, a few extreme branches (misaligned sequences, long-branch artefacts, wrong-topology neighbours) dominate. Use the MAD-based value as the baseline and SD as an upper bound.
3. **Check heteroscedasticity.** If `σ` is clearly larger in the shortest-branch bin, a single constant σ overstates noise for long branches; consider σ as a function of `L_ref` (binned or `σ = a + b/√L_ref`).
4. **Terminal vs internal.** Expect `σ_internal > σ_external` (short internal branches are poorly estimated). Keep them separate.
5. **Species-level structure.** If the species-demeaned `σ_external` is close to `σ_external`, pooled σ is adequate. If species differ a lot (`σ_s` spread, large `|mean_E|`), use `--mode species`, and consider adding a per-species offset.
6. **Do not double count estimation error.** `E` contains true rate variation *and* gene-tree estimation error. If simulated trees are fed into sequence simulation and re-inferred, inference adds noise a second time. Treat the empirical σ as an **upper bound** and use the self-consistency check (Step 4 in `run_all.sh`), extended to re-infer trees from simulated alignments of matching length, and shrink σ until the re-estimated `σ_external`, `σ_internal` match the empirical values.
7. Final recommendation = baseline (MAD- or SD-based, after 2–6) plus a sensitivity grid `--scale 0.5, 1, 1.5`.

## 6. Validation and sanity checks

- **Self-consistency:** rerun steps 1–2 on simulated trees; recovered `σ_external`, `σ_internal`, `σ_gene` should match the inputs. (On a synthetic test with known σ = 0.30 / 0.50 / 0.20, the pipeline recovered 0.335 / 0.49 / 0.206 after 20% of genes were given deliberate tip swaps; the MAD estimate for external was 0.29, illustrating why the robust check matters.)
- **Split recovery:** distribution of `frac_ref_splits_recovered`. Low recovery means few internal branches and fragile `σ_internal`.
- **Censoring rate:** fraction of branches below `--min-length`. If large, the log-scale SD is biased low; re-run with a different floor and compare.
- **`M_g` method:** compare median vs total-length versions.
- **Visual:** `diagnostics.png` (histograms of `E`, `E` vs reference length, per-species `σ_s`).

## 7. Caveats

- `E` is on the log scale, so `σ` is multiplicative variation. With `L = L_ref·exp(M+ε)` the *median* ratio is preserved; use `--mean-preserving` if the *mean* must be preserved.
- Because `M_g` is estimated from the same branches, residual SD is slightly underestimated for small genes (≈ 60–100 branches here, so minor).
- Reference terminal lengths of species whose sister was pruned include a merged internal length. This is correct for the pruned tree but makes `L_ref` gene-dependent.
- Branch estimates from the same gene are correlated; use gene-level bootstraps, not branch-level.
- Taxa with a bad reference length (e.g. very short) have inflated log ratios; inspect species with extreme `σ_s` or `mean_E` before using them.

## 8. Deliverables

`species_terminal_stats.csv`, `pooled_sigma.csv`, `sigma_by_length_bin.csv`, `sigma_by_clade_size.csv`, `gene_summary.csv`, `sim_params.json`, `diagnostics.png`, 100 simulated trees (+ `truth.csv`) per scenario, and the self-consistency check tables.
