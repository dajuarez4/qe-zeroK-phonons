#!/usr/bin/env python3
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
