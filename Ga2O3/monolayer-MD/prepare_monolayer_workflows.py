#!/usr/bin/env python3
"""Build 5x3x1 QE relaxation and 300 K MD workflows from monolayer FDF files."""

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent
REPEATS = (5, 3, 1)
VACUUM_C = 30.0


def read_fdf(path: Path):
    text = path.read_text()

    def block(name):
        match = re.search(
            rf"%block\s+{re.escape(name)}\s*\n(.*?)%endblock\s+{re.escape(name)}",
            text,
            re.I | re.S,
        )
        if not match:
            raise ValueError(f"Missing {name} block in {path}")
        return [line.split() for line in match.group(1).splitlines() if line.strip()]

    cell = [[float(x) for x in row[:3]] for row in block("LatticeVectors")]
    atoms = []
    for row in block("AtomicCoordinatesAndAtomicSpecies"):
        species = "O" if int(row[3]) == 1 else "Ga"
        atoms.append((species, *(float(x) for x in row[:3])))
    return cell, atoms


def supercell(cell, atoms):
    na, nb, nc = REPEATS
    old_c = cell[2][2]
    z_center = 0.5 * (min(atom[3] for atom in atoms) + max(atom[3] for atom in atoms))
    out_cell = [
        [na * value for value in cell[0]],
        [nb * value for value in cell[1]],
        [0.0, 0.0, VACUUM_C],
    ]
    out_atoms = []
    for ia in range(na):
        for ib in range(nb):
            for ic in range(nc):
                for symbol, x, y, z in atoms:
                    # Preserve Cartesian layer thickness when increasing the vacuum.
                    z_30 = 0.5 + (z - z_center) * old_c / VACUUM_C
                    out_atoms.append(
                        (symbol, (x + ia) / na, (y + ib) / nb, (z_30 + ic) / nc)
                    )
    return out_cell, out_atoms


def positions_text(atoms):
    return "\n".join(
        f"{symbol:2s}  {x: .12f}  {y: .12f}  {z: .12f}"
        for symbol, x, y, z in atoms
    )


def cell_text(cell):
    return "\n".join("  " + "  ".join(f"{x: .12f}" for x in row) for row in cell)


def qe_input(kind, label, cell, atoms):
    if kind == "relax":
        control = """  calculation   = 'relax'
  prefix        = '{label}-relax'
  outdir        = './tmp'
  pseudo_dir    = '../../pseudo'
  verbosity     = 'high'
  restart_mode  = 'from_scratch'
  nstep         = 150
  tstress       = .true.
  tprnfor       = .true.
  etot_conv_thr = 1.0d-5
  forc_conv_thr = 1.0d-4""".format(label=label)
        electrons_conv = "1.0d-10"
        ions = """&IONS
  ion_dynamics     = 'bfgs'
  trust_radius_ini = 0.05
  trust_radius_max = 0.20
  trust_radius_min = 1.0d-4
/"""
    else:
        control = """  calculation  = 'md'
  prefix       = '{label}-300K'
  outdir       = './tmp'
  pseudo_dir   = '../../pseudo'
  verbosity    = 'high'
  restart_mode = 'from_scratch'
  nstep        = 10000
  dt           = 20.0
  tstress      = .true.
  tprnfor      = .true.
  disk_io      = 'low'""".format(label=label)
        electrons_conv = "1.0d-8"
        ions = """&IONS
  ion_dynamics      = 'verlet'
  ion_temperature   = 'svr'
  tempw             = 300.0
  nraise            = 200
  pot_extrapolation = 'second_order'
  wfc_extrapolation = 'second_order'
/"""

    return f"""&CONTROL
{control}
/

&SYSTEM
  ibrav           = 0
  nat             = {len(atoms)}
  ntyp            = 2
  input_dft       = 'PBE'
  ecutwfc         = 80
  ecutrho         = 640
  occupations     = 'fixed'
  assume_isolated = '2D'
  nosym           = .true.
/

&ELECTRONS
  electron_maxstep = 200
  conv_thr         = {electrons_conv}
  mixing_beta      = 0.2
  diagonalization  = 'david'
/

{ions}

ATOMIC_SPECIES
Ga  69.723  Ga.pbe-dn-kjpaw_psl.1.0.0.UPF
O   15.999  O.pbe-n-kjpaw_psl.1.0.0.UPF

CELL_PARAMETERS angstrom
{cell_text(cell)}

ATOMIC_POSITIONS crystal
{positions_text(atoms)}

K_POINTS automatic
2 2 1 0 0 0
"""


RUN_TEMPLATE = """#!/bin/bash
#SBATCH -J {job}
#SBATCH -A DMR26002
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 4
#SBATCH -t 48:00:00
#SBATCH -o %x.%j.out
#SBATCH -e %x.%j.err

set -eo pipefail
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

cd "$SLURM_SUBMIT_DIR"
mkdir -p tmp
module purge
module load intel/19.1.1 impi/19.0.9

# QE 7.5 installation already verified in the user's Lonestar6 environment.
QE_BIN="${{QE_BIN:-${{SCRATCH}}/apps/qe-7.5-atomic-fixed-ntyp128/bin/pw.x}}"
if [[ ! -x "$QE_BIN" ]]; then
  echo "ERROR: QE executable not found or not executable: $QE_BIN" >&2
  exit 1
fi
{pre}
ibrun "$QE_BIN" -in {input_name} > {output_name}
grep -q "JOB DONE." {output_name}
{post}
"""


EXTRACTOR = r'''#!/usr/bin/env python3
"""Replace the MD template geometry with the final converged QE relax geometry."""
from pathlib import Path
import re

NAT = 75
HERE = Path(__file__).resolve().parent
RELAX_OUT = HERE.parent / "01-relax" / "relax.out"
TEMPLATE = HERE / "md.template.in"
TARGET = HERE / "md.in"

text = RELAX_OUT.read_text(errors="replace")
if "JOB DONE." not in text or "End of BFGS Geometry Optimization" not in text:
    raise SystemExit("Relaxation is absent or did not converge.")
if "maximum number of steps has been reached" in text:
    raise SystemExit("Relaxation reached nstep without force convergence.")

matches = re.findall(
    r"ATOMIC_POSITIONS\s*\(([^)]+)\)\s*\n"
    r"((?:\s*[A-Za-z][A-Za-z0-9]*\s+[-+0-9.EeDd]+\s+[-+0-9.EeDd]+\s+[-+0-9.EeDd]+[^\n]*\n){75})",
    text,
)
if not matches:
    raise SystemExit("Could not find the final 75-atom ATOMIC_POSITIONS block.")
unit, block = matches[-1]
if "crystal" not in unit.lower():
    raise SystemExit(f"Expected crystal coordinates, found {unit}.")

clean = []
for line in block.splitlines()[:NAT]:
    fields = line.split()
    clean.append(
        f"{fields[0]:2s}  {float(fields[1].replace('D','E')): .12f}  "
        f"{float(fields[2].replace('D','E')): .12f}  {float(fields[3].replace('D','E')): .12f}"
    )

template = TEMPLATE.read_text()
updated, count = re.subn(
    r"(ATOMIC_POSITIONS\s+crystal\s*\n).*?(\n\s*K_POINTS)",
    lambda m: m.group(1) + "\n".join(clean) + m.group(2),
    template,
    count=1,
    flags=re.S,
)
if count != 1:
    raise SystemExit("Could not replace the MD template coordinates.")
TARGET.write_text(updated)
print(f"Wrote {TARGET} from the converged relaxation.")
'''


def build(name, fdf):
    cell, atoms = read_fdf(fdf)
    cell, atoms = supercell(cell, atoms)
    if len(atoms) != 75:
        raise ValueError(f"Expected 75 atoms for {name}, found {len(atoms)}")

    top = ROOT / f"Ga2O3-ML-{name}-5x3x1-300K"
    relax_dir = top / "01-relax"
    md_dir = top / "02-md-300K"
    relax_dir.mkdir(parents=True, exist_ok=True)
    md_dir.mkdir(parents=True, exist_ok=True)

    label = f"Ga2O3-ML-{name}-5x3x1"
    (relax_dir / "relax.in").write_text(qe_input("relax", label, cell, atoms))
    smoke = qe_input("relax", f"{label}-smoke", cell, atoms)
    smoke = smoke.replace("calculation   = 'relax'", "calculation   = 'scf'")
    smoke = re.sub(r"^\s*(nstep|etot_conv_thr|forc_conv_thr)\s*=.*\n", "", smoke, flags=re.M)
    smoke = re.sub(r"\n&IONS\n.*?\n/\n", "", smoke, flags=re.S)
    smoke = smoke.replace("2 2 1 0 0 0", "1 1 1 0 0 0")
    (relax_dir / "smoke_test.in").write_text(smoke)
    (md_dir / "md.template.in").write_text(qe_input("md", label, cell, atoms))
    (md_dir / "prepare_md.py").write_text(EXTRACTOR)
    (relax_dir / "run_relax.sbatch").write_text(
        RUN_TEMPLATE.format(
            job=f"ML{name}-relax", pre="", input_name="relax.in", output_name="relax.out",
            post='grep -q "End of BFGS Geometry Optimization" relax.out',
        )
    )
    smoke_job = f"""#!/bin/bash
#SBATCH -J ML{name}-smoke
#SBATCH -A DMR26002
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 4
#SBATCH -t 00:15:00
#SBATCH -o smoke-test.%j.out
#SBATCH -e smoke-test.%j.err

set -eo pipefail
export OMP_NUM_THREADS=4
export MKL_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

cd "$SLURM_SUBMIT_DIR"
mkdir -p tmp-smoke
module purge
module load intel/19.1.1 impi/19.0.9

QE_BIN="${{QE_BIN:-${{SCRATCH}}/apps/qe-7.5-atomic-fixed-ntyp128/bin/pw.x}}"
if [[ ! -x "$QE_BIN" ]]; then
  echo "ERROR: QE executable not found or not executable: $QE_BIN" >&2
  exit 1
fi

echo "QE executable: $QE_BIN"
ibrun "$QE_BIN" -in smoke_test.in > smoke_test.qe.out
grep -q "JOB DONE." smoke_test.qe.out
grep -q "convergence has been achieved" smoke_test.qe.out
echo "PASS: 75-atom {name} SCF smoke test completed successfully."
"""
    # Keep smoke-test scratch separate from the production relaxation.
    smoke = (relax_dir / "smoke_test.in").read_text().replace("outdir        = './tmp'", "outdir        = './tmp-smoke'")
    (relax_dir / "smoke_test.in").write_text(smoke)
    (relax_dir / "run_smoke_test.sbatch").write_text(smoke_job)
    (md_dir / "run_md.sbatch").write_text(
        RUN_TEMPLATE.format(
            job=f"ML{name}-MD", pre="python3 prepare_md.py", input_name="md.in",
            output_name="md.out", post="",
        )
    )
    submit = """#!/bin/bash
set -euo pipefail
root="$(cd "$(dirname "$0")" && pwd)"
relax_submit=$(cd "$root/01-relax" && sbatch --parsable run_relax.sbatch)
# On federated Slurm systems --parsable may return JOBID;CLUSTER.
relax_id="${relax_submit%%;*}"
if [[ ! "$relax_id" =~ ^[0-9]+$ ]]; then
    echo "ERROR: could not extract a numeric relaxation job ID from: $relax_submit" >&2
    exit 1
fi
md_submit=$(cd "$root/02-md-300K" && sbatch --parsable --dependency="afterok:${relax_id}" run_md.sbatch)
md_id="${md_submit%%;*}"
echo "relax job: $relax_id"
echo "MD job:    $md_id (afterok:$relax_id)"
"""
    (top / "submit_workflow.sh").write_text(submit)


build("001", ROOT / "fdf_stru" / "Ga2O3-GGA-Bands-ml-001.fdf")
build("100", ROOT / "fdf_stru" / "Ga2O3-GGA-Bands-ml-100.fdf")
print("Prepared ML-001 and ML-100 5x3x1 (75-atom) workflows.")
