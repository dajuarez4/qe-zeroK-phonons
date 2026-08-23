#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
outputs=(displacements/disp-*/alpha_Fe2O3.fd.scf.out)
if [ "${#outputs[@]}" -ne 6 ]; then
  echo "Expected six QE outputs, found ${#outputs[@]}."; exit 1
fi
for output in "${outputs[@]}"; do
  grep -q 'convergence has been achieved' "$output" || { echo "Unconverged: $output"; exit 1; }
  grep -q 'Forces acting on atoms' "$output" || { echo "No forces: $output"; exit 1; }
  grep -q 'JOB DONE' "$output" || { echo "Incomplete: $output"; exit 1; }
done
python3 ../check_force_outputs.py
phonopy --qe -f "${outputs[@]}"
phonopy-load phonopy_disp.yaml --writefc
phonopy-load phonopy_disp.yaml --config band.conf --band-format hdf5
phonopy-load phonopy_disp.yaml --mesh 20 20 20 --gc --dos
echo "Force constants, bands, and DOS completed."
