#!/usr/bin/env python3
"""Assemble tabGAP forces and calculate 0 K Phonopy dispersion and DOS."""

from pathlib import Path
import numpy as np
import phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections

HERE = Path(__file__).resolve().parent


def read_forces(path):
    lines = path.read_text().splitlines()
    start = next(i for i,x in enumerate(lines) if x.startswith("ITEM: ATOMS"))
    cols = lines[start].split()[2:]
    rows = [dict(zip(cols,x.split())) for x in lines[start+1:]]
    rows.sort(key=lambda r:int(r["id"]))
    return [[float(r[c]) for c in ("fx","fy","fz")] for r in rows]


def main():
    ph = phonopy.load(HERE/"phonopy_disp.yaml")
    files = sorted((HERE/"forces").glob("force-*.dump"))
    expected = len(ph.supercells_with_displacements)
    if len(files) != expected:
        raise RuntimeError(f"Expected {expected} force files, found {len(files)}")
    ph.forces = np.array([read_forces(p) for p in files])
    ph.produce_force_constants(fc_calculator="symfc")
    ph.symmetrize_force_constants()
    ph.save(HERE/"phonopy_0K_tabgap.yaml", settings={"force_constants": True})

    points = [[[0,0,0],[0.5,0,0],[0.5,0.5,0],[0,0.5,0],[0,0,0]]]
    bands, connections = get_band_qpoints_and_path_connections(points, npoints=101)
    ph.run_band_structure(bands, path_connections=connections,
                          labels=[r"$\\Gamma$","X","S","Y",r"$\\Gamma$"])
    bd = ph.get_band_structure_dict()
    xparts, fparts = bd["distances"], bd["frequencies"]
    x0=0.0; output=[]
    for x,f in zip(xparts,fparts):
        xx=x-x[0]+x0; x0=xx[-1]
        output.append(np.column_stack([xx,f]))
    arr=np.vstack(output)
    np.savetxt(HERE/"dispersion_0K_tabgap.dat",arr)

    ph.run_mesh([24,24,1], with_eigenvectors=False, is_mesh_symmetry=True)
    ph.run_total_dos()
    dos=ph.get_total_dos_dict()
    np.savetxt(HERE/"dos_0K_tabgap.dat",np.column_stack([dos["frequency_points"],dos["total_dos"]]))
    f=np.vstack(fparts)
    print(f"Finite-displacement minimum: {f.min():.9f} THz")
    print(f"Negative values: {(f<0).sum()} across {(f<0).any(axis=1).sum()} of {len(f)} q points")


if __name__ == "__main__": main()
