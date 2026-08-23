#!/usr/bin/env python3
"""Prepare matched AFM/FM finite-displacement phonons at the BM equilibrium."""

from __future__ import annotations

from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
import yaml


ROOT = Path(__file__).resolve().parent
TARGET = ROOT / "spin_phonon_BM_equilibrium"
V0 = 106.261887  # A^3 per 10-atom primitive cell, nine-point BM3 fit
V104 = 104.6499400059482
V106 = 106.66245139092278
AFM_MOMENTS = (0, 0, 0, 0, 0, 0, 4, 4, -4, -4)


def read_structure(label: str) -> tuple[np.ndarray, list[tuple[str, np.ndarray]]]:
    text = (ROOT / f"phonopy_{label}_2x2x2" / "unitcell.in").read_text()
    cell_match = re.search(r"CELL_PARAMETERS angstrom\n((?:[^\n]+\n){3})", text)
    pos_match = re.search(r"ATOMIC_POSITIONS crystal\n((?:[^\n]+\n){10})", text)
    if not cell_match or not pos_match:
        raise RuntimeError(f"Could not parse the {label} unit cell")
    cell = np.array(
        [[float(x) for x in row.split()] for row in cell_match.group(1).splitlines()]
    )
    positions = []
    for row in pos_match.group(1).splitlines():
        fields = row.split()
        positions.append((fields[0], np.array([float(x) for x in fields[1:4]])))
    return cell, positions


def equilibrium_structure() -> tuple[np.ndarray, list[tuple[str, np.ndarray]]]:
    cell104, pos104 = read_structure("V104")
    _, pos106 = read_structure("V106")
    fraction = (V0 - V104) / (V106 - V104)
    cell = cell104 * (V0 / abs(np.linalg.det(cell104))) ** (1.0 / 3.0)
    positions = []
    for (name104, xyz104), (name106, xyz106) in zip(pos104, pos106):
        if name104 != name106:
            raise RuntimeError("V104 and V106 atom ordering differs")
        delta = xyz106 - xyz104
        delta -= np.rint(delta)
        positions.append((name104, (xyz104 + fraction * delta) % 1.0))
    if not np.isclose(abs(np.linalg.det(cell)), V0, atol=1e-6):
        raise RuntimeError("Failed to construct the BM equilibrium volume")
    return cell, positions


def qe_structure(cell: np.ndarray, positions: list[tuple[str, np.ndarray]]) -> str:
    cell_rows = "\n".join("  " + "  ".join(f"{x: .12f}" for x in row) for row in cell)
    atom_rows = "\n".join(
        f"{name:<4s}  {xyz[0]:.10f}  {xyz[1]:.10f}  {xyz[2]:.10f}"
        for name, xyz in positions
    )
    return f"""ATOMIC_SPECIES
O       15.999  O.pbe-n-kjpaw_psl.1.0.0.UPF
Fe1     55.845  Fe.pbe-spn-kjpaw_psl.1.0.0.UPF
Fe2     55.845  Fe.pbe-spn-kjpaw_psl.1.0.0.UPF

CELL_PARAMETERS angstrom
{cell_rows}

ATOMIC_POSITIONS crystal
{atom_rows}
"""


def unitcell_input(structure: str, state: str) -> str:
    total = "0.0" if state == "AFM" else "20.0"
    second = "-0.80" if state == "AFM" else " 0.80"
    return f"""&CONTROL
  calculation = 'scf'
  prefix      = 'alpha-Fe2O3-BMeq-{state}'
  outdir      = './tmp'
  pseudo_dir  = '../../pseudo'
  verbosity   = 'high'
  tstress     = .true.
  tprnfor     = .true.
/

&SYSTEM
  ibrav             = 0
  nat               = 10
  ntyp              = 3
  input_dft         = 'PBE'
  ecutwfc           = 80
  ecutrho           = 640
  occupations       = 'fixed'
  nspin             = 2
  tot_magnetization = {total}
  starting_magnetization(2) =  0.80
  starting_magnetization(3) = {second}
  lda_plus_u        = .true.
  Hubbard_U(2)      = 4.0
  Hubbard_U(3)      = 4.0
  U_projection_type = 'ortho-atomic'
/

&ELECTRONS
  electron_maxstep = 300
  conv_thr         = 1.0d-12
  mixing_beta      = 0.15
  mixing_mode      = 'local-TF'
  diagonalization  = 'david'
/

{structure}
"""


def force_input_header(state: str) -> str:
    total = "0.0" if state == "AFM" else "160.0"
    second = "-0.80" if state == "AFM" else " 0.80"
    return f"""&CONTROL
  calculation = 'scf'
  prefix      = 'alpha-Fe2O3-BMeq-{state}-fd'
  outdir      = './tmp'
  pseudo_dir  = '../../../../pseudo'
  verbosity   = 'high'
  tstress     = .true.
  tprnfor     = .true.
  disk_io     = 'low'
/

&SYSTEM
  ibrav             = 0
  nat               = 80
  ntyp              = 3
  input_dft         = 'PBE'
  ecutwfc           = 80
  ecutrho           = 640
  occupations       = 'fixed'
  nspin             = 2
  tot_magnetization = {total}
  starting_magnetization(2) =  0.80
  starting_magnetization(3) = {second}
  lda_plus_u        = .true.
  Hubbard_U(2)      = 4.0
  Hubbard_U(3)      = 4.0
  U_projection_type = 'ortho-atomic'
  nosym             = .true.
  noinv             = .true.
/

&ELECTRONS
  electron_maxstep = 400
  conv_thr         = 1.0d-12
  mixing_beta      = 0.10
  mixing_mode      = 'local-TF'
  mixing_ndim      = 12
  diagonalization  = 'david'
/
"""


def prepare_script(state: str) -> str:
    header = repr(force_input_header(state))
    return f'''#!/usr/bin/env python3
"""Build {state} QE force inputs from the common displaced supercells."""

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent
HEADER = {header}
K_POINTS = "\\nK_POINTS automatic\\n4 4 4 0 0 0\\n"

sources = sorted(ROOT.glob("supercell-[0-9][0-9][0-9].in"))
if len(sources) != 6:
    raise SystemExit(f"Expected six displaced supercells, found {{len(sources)}}")
for index, source in enumerate(sources, 1):
    source_text = source.read_text()
    if not re.search(r"nat\\s*=\\s*80", source_text):
        raise SystemExit(f"{{source}} is not an 80-atom supercell")
    structure = "\\n".join(
        line for line in source_text.splitlines()
        if not line.lstrip().startswith("!")
    ).strip()
    directory = ROOT / "displacements" / f"disp-{{index:03d}}"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "alpha_Fe2O3.fd.scf.in").write_text(
        HEADER + "\\n" + structure + K_POINTS
    )
print(f"Prepared {{len(sources)}} {state} displacement inputs")
'''


def batch_script(state: str) -> str:
    return f"""#!/bin/bash
#SBATCH --job-name=Fe2O3-BMeq-{state}-FD
#SBATCH --account=jakar_general
#SBATCH --qos=jakar_medium_general
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=38
#SBATCH --exclusive
#SBATCH --array=0-5%2
#SBATCH --output=%A_%a.log

set -euo pipefail
pwd
hostname
date
ulimit -s unlimited
set +u
source /opt/intel/oneapi/setvars.sh -ofi_internal=1 --force
set -u
module purge
module load gnu12
module load openmpi4/4.1.5
module load tbb
module load compiler-rt
module load mkl
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

DISPLACEMENTS=(disp-001 disp-002 disp-003 disp-004 disp-005 disp-006)
POINT=${{DISPLACEMENTS[$SLURM_ARRAY_TASK_ID]}}
cd "$SLURM_SUBMIT_DIR/displacements/$POINT"
mpirun --mca pml ob1 --mca btl ^openib,uct \\
  -np "$SLURM_NTASKS" /shared/quantum-espresso/bin/pw.x \\
  -npool 2 -ndiag 1 \\
  -in alpha_Fe2O3.fd.scf.in > alpha_Fe2O3.fd.scf.out
grep -q 'convergence has been achieved' alpha_Fe2O3.fd.scf.out
grep -q 'Forces acting on atoms' alpha_Fe2O3.fd.scf.out
grep -q 'JOB DONE' alpha_Fe2O3.fd.scf.out
date
"""


def collect_script() -> str:
    return """#!/bin/bash
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
"""


def force_output_checker() -> str:
    return '''#!/usr/bin/env python3
"""Validate the electronic and magnetic state of six displaced outputs."""

from pathlib import Path
import re

state = Path.cwd().name
targets = {"AFM": 0.0, "FM": 160.0}
if state not in targets:
    raise SystemExit("Run this checker from the AFM or FM directory")

outputs = sorted(Path("displacements").glob("disp-*/alpha_Fe2O3.fd.scf.out"))
if len(outputs) != 6:
    raise SystemExit(f"Expected six outputs, found {len(outputs)}")
for output in outputs:
    text = output.read_text(errors="replace")
    if "convergence has been achieved" not in text or "JOB DONE" not in text:
        raise SystemExit(f"Unconverged or incomplete: {output}")
    if "Forces acting on atoms" not in text:
        raise SystemExit(f"No final forces: {output}")
    moments = re.findall(r"total magnetization\\s+=\\s+([-0-9.]+)", text)
    if not moments:
        raise SystemExit(f"No total magnetization found: {output}")
    moment = float(moments[-1])
    if abs(moment - targets[state]) > 0.2:
        raise SystemExit(
            f"Wrong magnetic state in {output}: M={moment:.3f}, expected {targets[state]:.1f}"
        )
    print(f"{output.parent.name}: converged, M={moment:.3f} muB")
print(f"All six {state} force calculations passed.")
'''


def validation_batch() -> str:
    return """#!/bin/bash
#SBATCH --job-name=Fe2O3-BMeq-magcheck
#SBATCH --account=jakar_general
#SBATCH --qos=jakar_medium_general
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=38
#SBATCH --exclusive
#SBATCH --array=0-1
#SBATCH --output=%A_%a.log

set -euo pipefail
set +u
source /opt/intel/oneapi/setvars.sh -ofi_internal=1 --force
set -u
module purge
module load gnu12
module load openmpi4/4.1.5
module load tbb
module load compiler-rt
module load mkl
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
STATES=(AFM FM)
STATE=${STATES[$SLURM_ARRAY_TASK_ID]}
cd "$SLURM_SUBMIT_DIR"
mpirun --mca pml ob1 --mca btl ^openib,uct \\
  -np "$SLURM_NTASKS" /shared/quantum-espresso/bin/pw.x \\
  -npool 2 -ndiag 1 -in "$STATE.in" > "$STATE.out"
grep -q 'convergence has been achieved' "$STATE.out"
grep -q 'JOB DONE' "$STATE.out"
"""


def readme() -> str:
    return f"""# AFM-FM spin-phonon comparison at the BM equilibrium

This workflow compares collinear AFM (`++--`) and constrained FM (`++++`)
phonons at exactly the same structure.

- BM3 equilibrium volume: {V0:.6f} A^3 per 10-atom cell
- BM3 hexagonal parameters: a = 5.127825 A, c = 13.999151 A
- internal coordinates: periodic linear interpolation between the independently
  relaxed V104 and V106 structures
- model: PBE+U, U(Fe 3d) = 4 eV, ortho-atomic projectors
- supercells: 2x2x2 (80 atoms), common six central displacements
- supercell k mesh: 4x4x4

The BM point belongs to the existing fixed-shape EOS path. It is the energy-fit
minimum, not a new fully anisotropic zero-stress `vc-relax` structure.

## 1. Validate the two magnetic states

```bash
cd validation
sbatch run_validation.sbatch
python3 check_validation.py
```

Proceed only if both SCFs converged, AFM has total magnetization near 0, FM
has total magnetization near 20 muB per primitive cell, and the reported local
Fe moments preserve their intended signs.

## 2. Submit matched force calculations

```bash
cd ../AFM
sbatch run_displacements.sbatch
cd ../FM
sbatch run_displacements.sbatch
```

Each state has six calculations, limited to two simultaneous jobs. The FM
calculation is constrained to 160 muB in the 80-atom supercell.

## 3. Post-process after all jobs finish

```bash
cd AFM && bash collect_forces.sh
cd ../FM && bash collect_forces.sh
cd .. && python3 plot_afm_fm_comparison.py
```

The same displacement set and the AFM magnetic subgroup symmetry are used for
both fits. This deliberately avoids attributing different symmetry reduction
or different structures to spin-phonon coupling.
"""


def validation_checker() -> str:
    return '''#!/usr/bin/env python3
"""Check convergence, energies, and total moments of AFM and FM validation runs."""

from pathlib import Path
import re

for state, target in (("AFM", 0.0), ("FM", 20.0)):
    path = Path(f"{state}.out")
    text = path.read_text(errors="replace")
    if "convergence has been achieved" not in text or "JOB DONE" not in text:
        raise SystemExit(f"{state} did not finish with a converged SCF")
    energies = re.findall(r"!\\s+total energy\\s+=\\s+([-0-9.]+)\\s+Ry", text)
    moments = re.findall(r"total magnetization\\s+=\\s+([-0-9.]+)", text)
    if not energies or not moments:
        raise SystemExit(f"Could not parse {state} energy or magnetization")
    moment = float(moments[-1])
    if abs(moment - target) > 0.2:
        raise SystemExit(f"{state} total magnetization {moment:.3f} differs from {target:.1f}")
    print(f"{state}: E = {float(energies[-1]):.10f} Ry, M = {moment:.3f} muB")
print("Inspect the final site-resolved Fe moments in both outputs before submitting phonons.")
'''


def main() -> None:
    if TARGET.exists():
        raise FileExistsError(f"Refusing to overwrite {TARGET}")
    phonopy = shutil.which("phonopy")
    if not phonopy:
        raise SystemExit("phonopy was not found")

    cell, positions = equilibrium_structure()
    structure = qe_structure(cell, positions)
    TARGET.mkdir()
    (TARGET / "README.md").write_text(readme())
    (TARGET / "check_force_outputs.py").write_text(force_output_checker())

    validation = TARGET / "validation"
    validation.mkdir()
    for state in ("AFM", "FM"):
        text = unitcell_input(structure, state) + "\nK_POINTS automatic\n8 8 8 0 0 0\n"
        (validation / f"{state}.in").write_text(text)
    (validation / "run_validation.sbatch").write_text(validation_batch())
    (validation / "check_validation.py").write_text(validation_checker())

    template = ROOT / "phonopy_V104_2x2x2"
    afm = TARGET / "AFM"
    afm.mkdir()
    (afm / "unitcell.in").write_text(unitcell_input(structure, "AFM"))
    for name in ("band.conf", "mesh.conf"):
        shutil.copy2(template / name, afm / name)
    subprocess.run(
        [
            phonopy, "--qe", "-d", "--dim", "2", "2", "2",
            "--amplitude", "0.02", "--pm", "--magmom",
            *[str(x) for x in AFM_MOMENTS], "-c", "unitcell.in",
        ],
        cwd=afm,
        check=True,
    )
    generated = sorted(afm.glob("supercell-[0-9][0-9][0-9].in"))
    if len(generated) not in (6, 12):
        raise RuntimeError(f"Expected six or twelve generated displacements, found {len(generated)}")

    # Phonopy 2.22 interprets QE extended labels Fe1/Fe2 as He/Li and finds
    # only half of the magnetic operations recognized by the 2.28 workflow
    # used for the original volume series. Retain its correctly generated
    # equilibrium coordinates, but transplant them into the validated V104
    # 2.28 metadata and select the same six central displacements.
    current_yaml = yaml.safe_load((afm / "phonopy_disp.yaml").read_text())
    reference_yaml = yaml.safe_load((template / "phonopy_disp.yaml").read_text())
    for section in ("primitive_cell", "unit_cell", "supercell"):
        reference_yaml[section]["lattice"] = current_yaml[section]["lattice"]
        for reference_point, current_point in zip(
            reference_yaml[section]["points"], current_yaml[section]["points"]
        ):
            reference_point["coordinates"] = current_point["coordinates"]
    selected_indices = (0, 1, 2, 3, 8, 9) if len(generated) == 12 else tuple(range(6))
    reference_yaml["displacements"] = [
        current_yaml["displacements"][index] for index in selected_indices
    ]
    (afm / "phonopy_disp.yaml").write_text(
        yaml.safe_dump(reference_yaml, sort_keys=False, width=1000)
    )

    def restore_qe_labels(text: str) -> str:
        text = re.sub(r"(?m)^(\s*)He(\s+)", r"\1Fe1\2", text)
        return re.sub(r"(?m)^(\s*)Li(\s+)", r"\1Fe2\2", text)

    selected_text = [restore_qe_labels(generated[index].read_text()) for index in selected_indices]
    for path in generated:
        path.unlink()
    for index, text in enumerate(selected_text, 1):
        (afm / f"supercell-{index:03d}.in").write_text(text)
    (afm / "supercell.in").write_text(
        restore_qe_labels((afm / "supercell.in").read_text())
    )
    supercells = sorted(afm.glob("supercell-[0-9][0-9][0-9].in"))

    fm = TARGET / "FM"
    fm.mkdir()
    (fm / "unitcell.in").write_text(unitcell_input(structure, "FM"))
    for path in [afm / "phonopy_disp.yaml", afm / "supercell.in", *supercells]:
        shutil.copy2(path, fm / path.name)
    for name in ("band.conf", "mesh.conf"):
        shutil.copy2(template / name, fm / name)

    for state, directory in (("AFM", afm), ("FM", fm)):
        (directory / "prepare_qe_inputs.py").write_text(prepare_script(state))
        (directory / "run_displacements.sbatch").write_text(batch_script(state))
        (directory / "collect_forces.sh").write_text(collect_script())
        subprocess.run(["python3", "prepare_qe_inputs.py"], cwd=directory, check=True)

    for path in TARGET.rglob("*.sbatch"):
        path.chmod(0o755)
    for path in TARGET.rglob("*.sh"):
        path.chmod(0o755)
    for path in TARGET.rglob("*.py"):
        path.chmod(0o755)
    print(f"Prepared {TARGET.relative_to(ROOT)} at V={abs(np.linalg.det(cell)):.6f} A^3")


if __name__ == "__main__":
    main()
