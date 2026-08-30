#!/bin/bash

set -euo pipefail

PROJECT="/Users/dajuarez4/Documents/qe-zeroK-phonons/Ga2O3/bilayer-MD/Ga2O3-BL-001-4A-relax-4x2x1-300K"

cd "$PROJECT"

if [[ ! -s Ga2O3-BL-001-4x2x1.md.out ]]; then
    echo "ERROR: no se encontró Ga2O3-BL-001-4x2x1.md.out"
    exit 1
fi

TDEP_OUTDIR=TDEP_300K \
TDEP_CUTOFF=5.5 \
TDEP_BIN_DIR=/Users/dajuarez4/Documents/Fe/tdep/build/src \
PYTHON_BIN=/opt/anaconda3/bin/python \
bash run_tdep.sh

echo
echo "===== ÚLTIMAS LÍNEAS DEL AJUSTE TDEP ====="
tail -n 5 TDEP_300K/extract_forceconstants.log

echo
echo "===== ÚLTIMAS LÍNEAS DEL CÁLCULO DE FONONES ====="
tail -n 5 TDEP_300K/phonon_dispersion_relations.log

echo
echo "===== RESUMEN FINAL ====="
/opt/anaconda3/bin/python -m json.tool TDEP_300K/tdep_summary.json

echo
echo "Workflow TDEP completado."
echo "Gráfica: $PROJECT/TDEP_300K/Ga2O3-BL-001-4A-TDEP-diagnostic.png"
