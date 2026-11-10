#!/usr/bin/env bash
# Usage: bash run_all.sh reference.tre gene_trees_dir "*.treefile"
set -euo pipefail
REF=$1; GENES=$2; PAT=${3:-"*"}
S=$(dirname "$0")/scripts
pip install --quiet --break-system-packages dendropy pandas numpy matplotlib

# 1-2: empirical branch statistics and sigma estimates
python $S/01_extract_branch_stats.py --ref $REF --genes $GENES --pattern "$PAT" --out results/empirical
python $S/02_summarise_sigma.py --in-dir results/empirical --out results/empirical

# 3: simulate 100 trees (pooled sigma; then species-specific and robust variants)
P=results/empirical/sim_params.json; T=results/empirical/gene_taxa.tsv
python $S/03_simulate_trees.py --ref $REF --params $P --n 100 --taxa-sets $T --out results/sim_pooled
python $S/03_simulate_trees.py --ref $REF --params $P --n 100 --taxa-sets $T --mode species --out results/sim_species
python $S/03_simulate_trees.py --ref $REF --params $P --n 100 --taxa-sets $T --robust --out results/sim_robust

# 4: self-consistency check - re-estimate sigma from simulated trees; should recover the inputs
for d in sim_pooled sim_species sim_robust; do
  python $S/01_extract_branch_stats.py --ref $REF --genes results/$d --pattern "sim_*.tre" --out results/check_$d
  python $S/02_summarise_sigma.py --in-dir results/check_$d --out results/check_$d --n-boot 200
done
