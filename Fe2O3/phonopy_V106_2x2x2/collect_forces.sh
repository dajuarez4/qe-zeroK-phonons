#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")"

outputs=(displacements/disp-*/alpha_Fe2O3.fd.scf.out)
if [ "${#outputs[@]}" -ne 6 ]; then
  echo "Expected 6 QE outputs, but found ${#outputs[@]}."
  exit 1
fi

for output in "${outputs[@]}"; do
  if ! grep -q 'JOB DONE' "$output"; then
    echo "Incomplete calculation: $output"
    exit 1
  fi
  if ! grep -q 'convergence has been achieved' "$output"; then
    echo "Unconverged SCF calculation: $output"
    exit 1
  fi
done

phonopy --qe -f "${outputs[@]}"
phonopy-load phonopy_disp.yaml --writefc

echo "Created FORCE_SETS and force_constants.hdf5/ FORCE_CONSTANTS."
echo "Next: phonopy-load phonopy_disp.yaml --config band.conf --save"
