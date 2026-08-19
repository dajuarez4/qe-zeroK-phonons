#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/dajuarez4/bilayer_ga2o3
LAMMPS="$ROOT/bin/lmp-tabgap"
TDEP=/home/dajuarez4/tdep-o0/bin

for temperature in 60 70 80 90 100; do
    folder="$ROOT/tdep_${temperature}K_40x40"
    cd "$folder"
    echo "===== ${temperature} K: LAMMPS MD ====="
    if [[ ! -s "trajectory_${temperature}K.lammpstrj" || ! -s "final_${temperature}K.data" ]]; then
        /usr/bin/mpirun -np 6 "$LAMMPS" -in "in.md_${temperature}K" -log log.lammps
    else
        echo "Existing complete-looking MD outputs found; skipping MD."
    fi

    echo "===== ${temperature} K: prepare five TDEP frames ====="
    python3 "prepare_tdep_${temperature}K.py"

    echo "===== ${temperature} K: fit second-order TDEP ====="
    "$TDEP/extract_forceconstants" -rc2 5.5 -s 1 --firstorder > extract_forceconstants.log 2>&1
    cp outfile.forceconstant infile.forceconstant
    cp "$ROOT/tdep_10K_40x40/infile.qpoints_dispersion" .

    echo "===== ${temperature} K: dispersion and DOS ====="
    "$TDEP/phonon_dispersion_relations" --readpath --dos --qpoint_grid 24 24 1 > phonon_dispersion_relations.log 2>&1
    echo "===== ${temperature} K complete ====="
done
