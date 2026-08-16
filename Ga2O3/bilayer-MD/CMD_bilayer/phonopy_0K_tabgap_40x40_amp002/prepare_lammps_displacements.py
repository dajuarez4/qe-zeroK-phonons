#!/usr/bin/env python3
"""Convert Phonopy POSCAR displacement files to LAMMPS atomic data files."""

from pathlib import Path

HERE = Path(__file__).resolve().parent


def read_poscar(path):
    lines = path.read_text().splitlines()
    scale = float(lines[1])
    cell = [[scale*float(x) for x in lines[i].split()] for i in range(2, 5)]
    names = lines[5].split(); counts = [int(x) for x in lines[6].split()]
    mode_line = 7
    if lines[mode_line].strip().lower().startswith("s"):
        mode_line += 1
    direct = lines[mode_line].strip().lower().startswith("d")
    species = sum(([s]*n for s,n in zip(names, counts)), [])
    coords = [[float(x) for x in lines[i].split()[:3]] for i in range(mode_line+1, mode_line+1+len(species))]
    if direct:
        coords = [[sum(f[j]*cell[j][k] for j in range(3)) for k in range(3)] for f in coords]
    return cell, species, coords


def write_data(source, target):
    cell, species, xyz = read_poscar(source)
    a,b,c = cell
    assert abs(a[1])+abs(a[2])+abs(b[2])+abs(c[0])+abs(c[1]) < 1e-7
    with target.open("w") as out:
        out.write(f"Phonopy displacement from {source.name}\n\n{len(xyz)} atoms\n2 atom types\n\n")
        out.write(f"0.0 {a[0]:.16f} xlo xhi\n0.0 {b[1]:.16f} ylo yhi\n0.0 {c[2]:.16f} zlo zhi\n")
        out.write(f"{b[0]:.16f} 0.0 0.0 xy xz yz\n\nMasses\n\n1 69.723\n2 15.999\n\nAtoms # atomic\n\n")
        for i,(s,r) in enumerate(zip(species,xyz),1):
            out.write(f"{i} {1 if s=='Ga' else 2} {r[0]:.16f} {r[1]:.16f} {r[2]:.16f}\n")


def main():
    outdir = HERE/"lammps_data"; outdir.mkdir(exist_ok=True)
    files = sorted(HERE.glob("POSCAR-*"))
    for i,path in enumerate(files,1):
        write_data(path, outdir/f"disp-{i:03d}.data")
    print(f"Prepared {len(files)} LAMMPS displacement files")


if __name__ == "__main__": main()

