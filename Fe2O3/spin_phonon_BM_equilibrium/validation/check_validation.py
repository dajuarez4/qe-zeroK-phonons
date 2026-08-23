#!/usr/bin/env python3
"""Check convergence, energies, and total moments of AFM and FM validation runs."""

from pathlib import Path
import re

for state, target in (("AFM", 0.0), ("FM", 20.0)):
    path = Path(f"{state}.out")
    text = path.read_text(errors="replace")
    if "convergence has been achieved" not in text or "JOB DONE" not in text:
        raise SystemExit(f"{state} did not finish with a converged SCF")
    energies = re.findall(r"!\s+total energy\s+=\s+([-0-9.]+)\s+Ry", text)
    moments = re.findall(r"total magnetization\s+=\s+([-0-9.]+)", text)
    if not energies or not moments:
        raise SystemExit(f"Could not parse {state} energy or magnetization")
    moment = float(moments[-1])
    if abs(moment - target) > 0.2:
        raise SystemExit(f"{state} total magnetization {moment:.3f} differs from {target:.1f}")
    print(f"{state}: E = {float(energies[-1]):.10f} Ry, M = {moment:.3f} muB")
print("Inspect the final site-resolved Fe moments in both outputs before submitting phonons.")
