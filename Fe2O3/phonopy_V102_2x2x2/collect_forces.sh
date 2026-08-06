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
    # A transferred stdout may lack only QE's final timing/footer section.
    # Accept it when the SCF converged and the complete force block is present.
    if grep -q 'convergence has been achieved' "$output" &&
       grep -q 'Forces acting on atoms' "$output" &&
       grep -q 'Total force =' "$output"; then
      echo "Warning: accepting complete forces from truncated output: $output"
    else
      echo "Incomplete calculation: $output"
      exit 1
    fi
  fi
done

phonopy --qe -f "${outputs[@]}"
phonopy-load phonopy_disp.yaml --writefc

echo "Created FORCE_SETS and FORCE_CONSTANTS."
echo "Next: phonopy-load phonopy_disp.yaml --config band.conf --save"
