#!/usr/bin/env python3
"""Convert the five 16,000-atom 10 K snapshots to TDEP input."""

from pathlib import Path
import shutil
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
NAT=16000
KB=8.617333262145e-5


def read_poscar(path):
    s=path.read_text().splitlines(); scale=float(s[1])
    cell=np.array([[float(x) for x in s[i].split()] for i in range(2,5)])*scale
    counts=[int(x) for x in s[6].split()]; n=sum(counts)
    frac=np.array([[float(x) for x in s[i].split()[:3]] for i in range(8,8+n)])
    return cell,frac


def read_dump(path):
    with path.open() as inp:
        while True:
            line=inp.readline()
            if not line: return
            if line!="ITEM: TIMESTEP\n": continue
            step=int(inp.readline()); assert inp.readline().startswith("ITEM: NUMBER")
            assert int(inp.readline())==NAT; assert inp.readline().startswith("ITEM: BOX")
            for _ in range(3): inp.readline()
            cols=inp.readline().split()[2:]
            rows=[dict(zip(cols,inp.readline().split())) for _ in range(NAT)]
            yield step,rows


def cart_to_frac(x,cell):
    fz=x[2]/cell[2,2]
    fy=(x[1]-fz*cell[2,1])/cell[1,1]
    fx=(x[0]-fy*cell[1,0]-fz*cell[2,0])/cell[0,0]
    return np.array([fx,fy,fz])


def main():
    uc=ROOT/"lammps_6x4_cellrelax/infile.ucposcar"
    ss=ROOT/"phonopy_0K_tabgap_40x40_amp002/SPOSCAR"
    shutil.copy2(uc,HERE/"infile.ucposcar"); shutil.copy2(ss,HERE/"infile.ssposcar")
    cell,ideal=read_poscar(ss); nframes=0; temps=[]
    with (HERE/"infile.positions").open("w") as pos, \
         (HERE/"infile.forces").open("w") as frc, \
         (HERE/"infile.stat").open("w") as stat:
        for step,rows in read_dump(HERE/"trajectory_10K.lammpstrj"):
            rows.sort(key=lambda r:int(r["id"]))
            observed=np.array([cart_to_frac([float(r["xu"]),float(r["yu"]),float(r["zu"])],cell) for r in rows])
            delta=(observed-ideal+.5)%1-.5
            drift=delta.mean(axis=0); observed=(observed-drift)%1
            pe=sum(float(r["c_atom_pe"]) for r in rows)
            ke=sum(float(r["c_atom_ke"]) for r in rows)
            temp=2*ke/(3*(NAT-1)*KB); temps.append(temp)
            np.savetxt(pos,observed,fmt="%.15e")
            np.savetxt(frc,[[float(r[c]) for c in ("fx","fy","fz")] for r in rows],fmt="%.15e")
            stat.write(f"{nframes+1} {nframes*2500.0:.8f} {pe+ke:.12e} {pe:.12e} {ke:.12e} {temp:.8f} "+"0.0 "*7+"\n")
            nframes+=1
    (HERE/"infile.meta").write_text(f"{NAT} # N atoms\n{nframes} # N configurations\n2500.0 # sampling interval, fs\n10.0 # temperature, K\n")
    print(f"Prepared {nframes} frames; mean temperature {np.mean(temps):.4f} K")


if __name__=="__main__": main()
