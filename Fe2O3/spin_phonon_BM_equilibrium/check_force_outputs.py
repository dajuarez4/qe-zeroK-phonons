#!/usr/bin/env python3
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
    moments = re.findall(r"total magnetization\s+=\s+([-0-9.]+)", text)
    if not moments:
        raise SystemExit(f"No total magnetization found: {output}")
    moment = float(moments[-1])
    if abs(moment - targets[state]) > 0.2:
        raise SystemExit(
            f"Wrong magnetic state in {output}: M={moment:.3f}, "
            f"expected {targets[state]:.1f}"
        )
    print(f"{output.parent.name}: converged, M={moment:.3f} muB")
print(f"All six {state} force calculations passed.")
