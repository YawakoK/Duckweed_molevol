#!/bin/bash
# Orchestrates steps 02 -> 03 -> 05 -> 06 -> 07 of 02_phylogenomics.
# Step 01 (OrthoFinder) and 04 (species tree) are run separately.
#
#   OGDIR   OrthoFinder Orthogroups/ directory
#   CDSDIR  per-species CDS files (see 02_build_orthogroup_cds.py)
#   NPROC   parallel jobs for the per-OG steps
#
# Everything here is compute-heavy (3,615 gene trees, 3,615 aBSREL fits);
# the resulting per-OG alignments, trees and aBSREL JSONs are archived on
# Zenodo, and the summary tables they yield are shipped in data/ so that
# steps 08-09 (the figures) can be run without repeating this.
set -euo pipefail
OGDIR=${OGDIR:?set OGDIR to the OrthoFinder Orthogroups/ directory}
CDSDIR=${CDSDIR:?set CDSDIR to the per-species CDS directory}
export TREES=${TREES:-TreeEachOG}
export ABSREL=${ABSREL:-absrel/absrel_json}
WORK=${WORK:-absrel}
NPROC=${NPROC:-26}
mkdir -p "$WORK" "$TREES"

echo "############ 1. orthogroup set + per-OG CDS ############"
python3 02_build_orthogroup_cds.py --ogdir "$OGDIR" --cds "$CDSDIR" --out "$TREES" --work "$WORK" \
        | tee "$WORK/geneset.log"

echo "############ 2. codon alignment + ML gene tree ############"
ls -d "$TREES"/OG* | xargs -n1 basename | sort > "$WORK/all_OG_list.txt"
xargs -a "$WORK/all_OG_list.txt" -n1 -P "$NPROC" bash 03_align_and_tree_one_og.sh \
      > "$WORK/build_progress.log" 2>&1 || true

echo "############ 3. aBSREL ############"
xargs -a "$WORK/all_OG_list.txt" -n1 -P "$NPROC" bash 05_absrel_one_og.sh \
      > "$WORK/absrel_progress.log" 2>&1 || true

echo "############ 4. summaries ############"
python3 06_summarise_absrel.py --json "$ABSREL" --trees "$TREES" --out output | tee output/absrel_summary.log
python3 07_roottotip_dnds.py   --json "$ABSREL" --per-branch output/absrel_per_branch.csv --out output
echo "############ PIPELINE COMPLETE $(date '+%F %H:%M') ############"
