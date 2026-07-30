#!/bin/bash

set -eo pipefail

TDEP_BIN_DIR="${TDEP_BIN_DIR:-/Users/dajuarez4/Documents/Fe/tdep/build/src}"

python3 qe_supercell_to_tdep.py \
  --unit-input ../Ga2O3-BL-001-300K/Ga2O3-BL-001.md.in \
  --supercell-input Ga2O3-BL-001-4x2x1.md.in \
  --qe-output Ga2O3-BL-001-4x2x1.md.out \
  --outdir TDEP_300K

cd TDEP_300K

"$TDEP_BIN_DIR/extract_forceconstants/extract_forceconstants" \
  --secondorder_cutoff 2.5 \
  --temperature 300 \
  --firstorder \
  --norotational \
  --nohuang > extract_forceconstants.log 2>&1

cp outfile.forceconstant infile.forceconstant

"$TDEP_BIN_DIR/phonon_dispersion_relations/phonon_dispersion_relations" \
  --readpath --unit thz > phonon_dispersion_relations.log 2>&1
