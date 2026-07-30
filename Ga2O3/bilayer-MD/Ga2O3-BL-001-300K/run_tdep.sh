#!/bin/bash

set -eo pipefail

TDEP_BIN_DIR="${TDEP_BIN_DIR:-/Users/dajuarez4/Documents/Fe/tdep/build/src}"
WORK_DIR="${1:-TDEP_300K}"

python3 qe_md_to_tdep.py \
  Ga2O3-BL-001.md.in Ga2O3-BL-001.md.out -o "$WORK_DIR"

cd "$WORK_DIR"

"$TDEP_BIN_DIR/extract_forceconstants/extract_forceconstants" \
  --secondorder_cutoff 2.5 \
  --temperature 300 \
  --firstorder \
  --norotational \
  --nohuang > extract_forceconstants.log 2>&1

cp outfile.forceconstant infile.forceconstant

"$TDEP_BIN_DIR/phonon_dispersion_relations/phonon_dispersion_relations" \
  --readpath \
  --unit thz > phonon_dispersion_relations.log 2>&1
