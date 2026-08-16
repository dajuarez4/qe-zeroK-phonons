#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
a=np.loadtxt(HERE/"dispersion_0K_tabgap.dat")
d=np.loadtxt(HERE/"dos_0K_tabgap.dat")
x=a[:,0]; f=a[:,1:]
bound=[x[0],x[100],x[201],x[302],x[403]]

fig,(ax,ad)=plt.subplots(1,2,figsize=(10.5,5.8),gridspec_kw={"width_ratios":[4,1]},sharey=True,constrained_layout=True)
for band in f.T: ax.plot(x,band,color="#7a4fa3",lw=0.85)
ax.axhline(0,color="black",lw=0.8)
for q in bound: ax.axvline(q,color="0.84",lw=0.7,zorder=0)
ax.set_xticks(bound,[r"$\Gamma$","X","S","Y",r"$\Gamma$"])
ax.set_xlim(bound[0],bound[-1]); ax.set_ylim(f.min()-0.4,f.max()+0.6)
ax.set_xlabel("Wave-vector path"); ax.set_ylabel("Frequency (THz)")
ax.grid(axis="y",color="0.93",lw=0.6)
ad.plot(d[:,1],d[:,0],color="#7a4fa3",lw=1.1)
ad.axhline(0,color="black",lw=0.8); ad.set_xlabel("DOS"); ad.set_xticks([])
fig.suptitle(r"Ga$_2$O$_3$ bilayer: 0 K tabGAP finite-displacement phonons\n6×4×1 supercell, 0.01 Å displacements",fontsize=14)
fig.savefig(HERE/"phonon_0K_tabgap_finite_displacement.png",dpi=220,bbox_inches="tight")
print(HERE/"phonon_0K_tabgap_finite_displacement.png")
