#!/bin/bash

set -eo pipefail

TDEP_BIN_DIR="${TDEP_BIN_DIR:-/Users/dajuarez4/Documents/Fe/tdep/build/src}"
PYTHON_BIN="${PYTHON_BIN:-/opt/anaconda3/bin/python}"
TDEP_OUTDIR="${TDEP_OUTDIR:-TDEP_300K}"
TDEP_CUTOFF="${TDEP_CUTOFF:-5.5}"

TDEP_OUTDIR="$TDEP_OUTDIR" "$PYTHON_BIN" prepare_tdep.py
cd "$TDEP_OUTDIR"

"$TDEP_BIN_DIR/extract_forceconstants/extract_forceconstants" \
  --secondorder_cutoff "$TDEP_CUTOFF" \
  --temperature 300 \
  --firstorder \
  --norotational \
  --nohuang > extract_forceconstants.log 2>&1

cp outfile.forceconstant infile.forceconstant

"$TDEP_BIN_DIR/phonon_dispersion_relations/phonon_dispersion_relations" \
  --readpath --unit thz > phonon_dispersion_relations.log 2>&1

cd ..
TDEP_OUTDIR="$TDEP_OUTDIR" TDEP_CUTOFF="$TDEP_CUTOFF" \
  "$PYTHON_BIN" plot_tdep_dispersion.py
