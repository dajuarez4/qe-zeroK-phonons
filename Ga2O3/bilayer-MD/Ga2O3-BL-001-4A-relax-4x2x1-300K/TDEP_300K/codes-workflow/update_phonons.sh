#!/usr/bin/env bash

# Refit the second-order force constants and update the TDEP phonons in the
# TDEP_300K directory that contains this codes-workflow folder.

set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
TDEP_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd -P)"
PROJECT_DIR="/Users/dajuarez4/Documents/qe-zeroK-phonons/Ga2O3/bilayer-MD/Ga2O3-BL-001-4A-relax-4x2x1-300K/"

TDEP_BIN_DIR="${TDEP_BIN_DIR:-/Users/dajuarez4/Documents/Fe/tdep/build/src}"
TDEP_CUTOFF="${TDEP_CUTOFF:-5.5}"
TDEP_TEMPERATURE="${TDEP_TEMPERATURE:-300}"
PYTHON_BIN="${PYTHON_BIN:-/opt/anaconda3/bin/python}"

EXTRACT_FC="${EXTRACT_FC:-$TDEP_BIN_DIR/extract_forceconstants/extract_forceconstants}"
PHONON_DISPERSION="${PHONON_DISPERSION:-$TDEP_BIN_DIR/phonon_dispersion_relations/phonon_dispersion_relations}"

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

check_file() {
    [[ -s "$TDEP_DIR/$1" ]] || die "missing or empty input: $TDEP_DIR/$1"
}

[[ "$TDEP_CUTOFF" =~ ^[0-9]+([.][0-9]+)?$ ]] || die "TDEP_CUTOFF must be a positive number"
[[ "$TDEP_TEMPERATURE" =~ ^[0-9]+([.][0-9]+)?$ ]] || die "TDEP_TEMPERATURE must be a positive number"
[[ -x "$EXTRACT_FC" ]] || die "extract_forceconstants is not executable: $EXTRACT_FC"
[[ -x "$PHONON_DISPERSION" ]] || die "phonon_dispersion_relations is not executable: $PHONON_DISPERSION"
[[ -x "$PYTHON_BIN" ]] || die "Python is not executable: $PYTHON_BIN"
[[ -f "$PROJECT_DIR/prepare_tdep.py" ]] || die "missing preparation code: $PROJECT_DIR/prepare_tdep.py"
check_file_from_project() {
    [[ -s "$PROJECT_DIR/$1" ]] || die "missing or empty source: $PROJECT_DIR/$1"
}
check_file_from_project "Ga2O3-BL-001-4x2x1.md.out"

RAW_MD_STEPS="$(awk '/Entering Dynamics:/ {n++} END {print n+0}' "$PROJECT_DIR/Ga2O3-BL-001-4x2x1.md.out")"
printf 'Raw QE MD steps found: %s\n' "$RAW_MD_STEPS"
printf 'Preparing fresh TDEP inputs from the current QE output...\n'
(
    cd "$PROJECT_DIR"
    TDEP_OUTDIR="$(basename -- "$TDEP_DIR")" "$PYTHON_BIN" prepare_tdep.py
)

for input in infile.ucposcar infile.ssposcar infile.positions infile.forces infile.stat infile.meta infile.qpoints_dispersion; do
    check_file "$input"
done

# TDEP's infile.meta begins with the supercell atom count followed by the
# number of configurations. Confirm that count against all trajectory files.
N_ATOMS="$(awk 'NF {print $1; exit}' "$TDEP_DIR/infile.meta")"
META_STEPS="$(awk 'NF {n++; if (n == 2) {print $1; exit}}' "$TDEP_DIR/infile.meta")"
[[ -n "$N_ATOMS" && -n "$META_STEPS" ]] || die "infile.meta does not contain atom and MD-step counts"
[[ "$N_ATOMS" =~ ^[1-9][0-9]*$ ]] || die "invalid atom count in infile.meta: $N_ATOMS"
[[ "$META_STEPS" =~ ^[1-9][0-9]*$ ]] || die "invalid MD-step count in infile.meta: $META_STEPS"

POSITION_LINES="$(awk 'NF {n++} END {print n+0}' "$TDEP_DIR/infile.positions")"
FORCE_LINES="$(awk 'NF {n++} END {print n+0}' "$TDEP_DIR/infile.forces")"
STAT_STEPS="$(awk 'NF {n++} END {print n+0}' "$TDEP_DIR/infile.stat")"
(( POSITION_LINES % N_ATOMS == 0 )) || die "infile.positions has an incomplete MD step"
(( FORCE_LINES % N_ATOMS == 0 )) || die "infile.forces has an incomplete MD step"
POSITION_STEPS=$((POSITION_LINES / N_ATOMS))
FORCE_STEPS=$((FORCE_LINES / N_ATOMS))

if (( META_STEPS != POSITION_STEPS || META_STEPS != FORCE_STEPS || META_STEPS != STAT_STEPS )); then
    die "MD-step counts disagree: meta=$META_STEPS, positions=$POSITION_STEPS, forces=$FORCE_STEPS, stat=$STAT_STEPS"
fi

printf 'TDEP directory: %s\n' "$TDEP_DIR"
printf 'Temperature:    %s K\n' "$TDEP_TEMPERATURE"
printf 'FC2 cutoff:    %s Angstrom\n' "$TDEP_CUTOFF"
printf 'Complete MD configurations prepared: %s (%s atoms each)\n' "$META_STEPS" "$N_ATOMS"

cd "$TDEP_DIR"

printf '\n[1/3] Fitting second-order force constants...\n'
"$EXTRACT_FC" \
    --secondorder_cutoff "$TDEP_CUTOFF" \
    --temperature "$TDEP_TEMPERATURE" \
    --firstorder \
    --norotational \
    --nohuang \
    > extract_forceconstants.log 2>&1

[[ -s outfile.forceconstant ]] || die "force-constant fit did not create outfile.forceconstant; inspect extract_forceconstants.log"
cp -f outfile.forceconstant infile.forceconstant

printf '[2/3] Calculating the phonon dispersion...\n'
"$PHONON_DISPERSION" \
    --readpath \
    --unit thz \
    > phonon_dispersion_relations.log 2>&1

[[ -s outfile.dispersion_relations ]] || die "phonon calculation did not create outfile.dispersion_relations; inspect phonon_dispersion_relations.log"

printf '[3/3] Updating the diagnostic plot...\n'
if [[ -f "$PROJECT_DIR/plot_tdep_dispersion.py" && -x "$PYTHON_BIN" ]]; then
    (
        cd "$PROJECT_DIR"
        TDEP_OUTDIR="$(basename -- "$TDEP_DIR")" \
        TDEP_CUTOFF="$TDEP_CUTOFF" \
        "$PYTHON_BIN" plot_tdep_dispersion.py
    )
else
    printf 'Skipped plot: plot_tdep_dispersion.py or Python was not found.\n'
fi

printf '\nDone. TDEP_300K phonons were updated successfully.\n'
printf 'Logs:\n  %s\n  %s\n' \
    "$TDEP_DIR/extract_forceconstants.log" \
    "$TDEP_DIR/phonon_dispersion_relations.log"
