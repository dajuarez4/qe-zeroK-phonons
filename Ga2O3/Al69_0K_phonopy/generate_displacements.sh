#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")"

if [ ! -f unitcell.in ]; then
  echo "unitcell.in is missing. Run extract_relaxed_structure.py after vc-relax."
  exit 1
fi

phonopy --qe -d --dim 1 2 1 --amplitude 0.02 --pm -c unitcell.in
python3 prepare_qe_inputs.py

count=$(find displacements -mindepth 1 -maxdepth 1 -type d -name 'disp-*' | wc -l | tr -d ' ')
if [ "$count" -lt 1 ]; then
  echo "No displacement directories were generated."
  exit 1
fi

echo "Generated $count central-difference force calculations."
echo "Submit with: sbatch --array=0-$((count - 1))%4 run_displacements.sbatch"
