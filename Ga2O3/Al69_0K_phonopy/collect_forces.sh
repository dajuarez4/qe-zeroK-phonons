#!/bin/bash

set -euo pipefail
cd "$(dirname "$0")"

mapfile -t outputs < <(find displacements -mindepth 2 -maxdepth 2 -type f -name 'Al69_Ga2O3.fd.scf.out' | sort)
expected=$(find displacements -mindepth 1 -maxdepth 1 -type d -name 'disp-*' | wc -l | tr -d ' ')
if [ "${#outputs[@]}" -ne "$expected" ]; then
  echo "Expected $expected QE outputs, found ${#outputs[@]}."
  exit 1
fi
for output in "${outputs[@]}"; do
  if ! grep -q 'JOB DONE' "$output"; then
    echo "Incomplete calculation: $output"
    exit 1
  fi
done

phonopy --qe -f "${outputs[@]}"
phonopy-load phonopy_disp.yaml --writefc
phonopy-load phonopy_disp.yaml --config band.conf --save
phonopy-load phonopy_disp.yaml --config mesh.conf
python3 plot_band_dos.py

echo "Created FORCE_SETS, FORCE_CONSTANTS, band, DOS, and thermal properties."
