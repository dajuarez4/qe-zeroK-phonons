#!/usr/bin/env python3
"""Project Phonopy FCs into hiPhive and enforce 2D rotational invariance."""

from pathlib import Path
import numpy as np
import phonopy
from ase.io import read
from hiphive import ClusterSpace, ForceConstantPotential, ForceConstants
from hiphive import enforce_rotational_sum_rules
from hiphive.utilities import extract_parameters
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections

HERE=Path(__file__).resolve().parent


def bands_for(ph):
    points=[[[0,0,0],[0.5,0,0],[0.5,0.5,0],[0,0.5,0],[0,0,0]]]
    bands,connections=get_band_qpoints_and_path_connections(points,npoints=101)
    ph.run_band_structure(bands,path_connections=connections)
    d=ph.get_band_structure_dict()
    return np.hstack(d["distances"]),np.vstack(d["frequencies"])


def read_forces(path):
    lines=path.read_text().splitlines()
    start=next(i for i,x in enumerate(lines) if x.startswith("ITEM: ATOMS"))
    cols=lines[start].split()[2:]
    rows=[dict(zip(cols,x.split())) for x in lines[start+1:]]
    rows.sort(key=lambda r:int(r["id"]))
    return [[float(r[c]) for c in ("fx","fy","fz")] for r in rows]

def main():
    ph=phonopy.load(HERE/"phonopy_disp.yaml")
    files=sorted((HERE/"forces").glob("force-*.dump"))
    ph.forces=np.array([read_forces(p) for p in files])
    ph.produce_force_constants(fc_calculator="symfc",calculate_full_force_constants=True)
    prim=read(HERE/"POSCAR")
    supercell=read(HERE/"SPOSCAR")
    raw=ForceConstants.from_arrays(supercell,ph.force_constants)

    # tabGAP's actual interaction cutoff is 7 A. A slightly smaller projection
    # cutoff avoids periodic image ambiguity at the boundary.
    cs=ClusterSpace(prim,[6.95])
    parameters=extract_parameters(raw,cs)
    fcp=ForceConstantPotential(cs,parameters)
    projected=fcp.get_force_constants(supercell)

    constrained=enforce_rotational_sum_rules(
        cs,parameters,["Huang","Born-Huang"],alpha=1e-8)
    fcp_rot=ForceConstantPotential(cs,constrained)
    rotational=fcp_rot.get_force_constants(supercell)
    fcp.write(str(HERE/"fcp_projected.fcp"))
    fcp_rot.write(str(HERE/"fcp_rotational.fcp"))

    results={}
    for name,fcs in (("raw",raw),("projected",projected),("rotational",rotational)):
        ph.force_constants=fcs.get_fc_array(order=2)
        x,f=bands_for(ph)
        np.savetxt(HERE/f"dispersion_{name}.dat",np.column_stack([x,f]))
        results[name]=(x,f)
        print(f"{name:10s}: min={f.min(): .9f} THz, negative={(f<0).sum()}")

    import matplotlib.pyplot as plt
    colors={"raw":"#555555","projected":"#d18f00","rotational":"#198754"}
    fig,axes=plt.subplots(1,2,figsize=(13,5.5),sharey=True,constrained_layout=True)
    for ax in axes:
        for name,(x,f) in results.items():
            for ib,b in enumerate(f.T):
                ax.plot(x,b,color=colors[name],lw=0.75,alpha=.75,
                        label=name if ib==0 else None)
        ax.axhline(0,color="black",lw=.8)
        ax.set_xlabel("Wave-vector path"); ax.grid(axis="y",color=".93",lw=.6)
    x=results["raw"][0]; bounds=[x[0],x[100],x[201],x[302],x[403]]
    for ax in axes:
        ax.set_xticks(bounds,[r"$\Gamma$","X","S","Y",r"$\Gamma$"])
        for q in bounds: ax.axvline(q,color=".85",lw=.6,zorder=0)
        ax.set_xlim(x[0],x[-1])
    axes[0].set_ylabel("Frequency (THz)"); axes[0].legend(frameon=False)
    axes[1].set_ylim(-1,5); axes[1].set_title("Acoustic-mode detail")
    axes[0].set_title("Full dispersion")
    fig.suptitle("0 K tabGAP phonons: rotational sum-rule convergence")
    fig.savefig(HERE/"phonon_rotational_sum_rules.png",dpi=220,bbox_inches="tight")


if __name__=="__main__": main()
