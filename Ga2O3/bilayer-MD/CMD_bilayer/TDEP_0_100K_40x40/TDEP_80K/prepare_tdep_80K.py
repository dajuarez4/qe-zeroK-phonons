#!/usr/bin/env python3
"""Convert five widely spaced frames from the 80 K trajectory to TDEP input."""

from pathlib import Path
import shutil
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
NAT = 16000
KB = 8.617333262145e-5
SELECT = {0, 6, 12, 19, 25}


def read_poscar(path):
    lines = path.read_text().splitlines()
    scale = float(lines[1])
    cell = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)]) * scale
    counts = [int(x) for x in lines[6].split()]
    n = sum(counts)
    frac = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(8, 8 + n)])
    return cell, frac


def read_dump(path):
    with path.open() as inp:
        iframe = 0
        while True:
            line = inp.readline()
            if not line:
                return
            if line != "ITEM: TIMESTEP\n":
                continue
            step = int(inp.readline())
            assert inp.readline().startswith("ITEM: NUMBER")
            assert int(inp.readline()) == NAT
            assert inp.readline().startswith("ITEM: BOX")
            for _ in range(3):
                inp.readline()
            columns = inp.readline().split()[2:]
            rows = [dict(zip(columns, inp.readline().split())) for _ in range(NAT)]
            yield iframe, step, rows
            iframe += 1


def cart_to_frac(x, cell):
    fz = x[2] / cell[2, 2]
    fy = (x[1] - fz * cell[2, 1]) / cell[1, 1]
    fx = (x[0] - fy * cell[1, 0] - fz * cell[2, 0]) / cell[0, 0]
    return np.array([fx, fy, fz])


uc = ROOT / "lammps_6x4_cellrelax/infile.ucposcar"
ss = ROOT / "phonopy_0K_tabgap_40x40_amp002/SPOSCAR"
shutil.copy2(uc, HERE / "infile.ucposcar")
shutil.copy2(ss, HERE / "infile.ssposcar")
cell, ideal = read_poscar(ss)
temperatures = []
nframes = 0

with (HERE / "infile.positions").open("w") as pos, \
     (HERE / "infile.forces").open("w") as frc, \
     (HERE / "infile.stat").open("w") as stat:
    for iframe, step, rows in read_dump(HERE / "trajectory_80K.lammpstrj"):
        if iframe not in SELECT:
            continue
        rows.sort(key=lambda row: int(row["id"]))
        observed = np.array([
            cart_to_frac([float(row["xu"]), float(row["yu"]), float(row["zu"])], cell)
            for row in rows
        ])
        delta = (observed - ideal + 0.5) % 1.0 - 0.5
        observed = (observed - delta.mean(axis=0)) % 1.0
        pe = sum(float(row["c_atom_pe"]) for row in rows)
        ke = sum(float(row["c_atom_ke"]) for row in rows)
        temperature = 2.0 * ke / (3.0 * (NAT - 1) * KB)
        temperatures.append(temperature)
        np.savetxt(pos, observed, fmt="%.15e")
        np.savetxt(frc, [[float(row[c]) for c in ("fx", "fy", "fz")] for row in rows], fmt="%.15e")
        stat.write(f"{nframes + 1} {step:.8f} {pe + ke:.12e} {pe:.12e} {ke:.12e} {temperature:.8f} " + "0.0 " * 7 + "\n")
        nframes += 1

(HERE / "infile.meta").write_text(
    f"{NAT} # N atoms\n{nframes} # N configurations\n"
    "1000.0 # trajectory sampling interval, fs\n80.0 # temperature, K\n"
)
print(f"Prepared {nframes} selected frames; mean temperature {np.mean(temperatures):.4f} K")
