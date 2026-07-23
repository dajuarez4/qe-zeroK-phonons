# Birch-Murnaghan equation of state for AFM alpha-Fe2O3

This is an independent equation-of-state workflow.  It does not edit or use
the completed outputs from the parent phonon calculation.

## What is calculated

Nine primitive-cell volumes are sampled symmetrically around the experimental
reference volume.  At each fixed volume, the ten atoms are relaxed while the
rhombohedral angle and hexagonal `c/a` ratio remain fixed.  This is the
standard isotropically scaled EOS path.

All points use exactly the same numerical and magnetic model:

- collinear `++--` antiferromagnetism
- `tot_magnetization=0`
- PBE+U with `U_eff=4.0 eV` on Fe 3d
- the same Fe and O PAW pseudopotentials as the phonon calculation
- 80/640 Ry cutoffs
- 8x8x8 k-point grid
- fixed occupations
- fixed cell volume with BFGS relaxation of internal coordinates

Do **not** change these inputs to `vc-relax`: every point would move toward the
same equilibrium volume and would no longer define an energy-volume curve.

## Volume points

The reference primitive volume is 100.624948 A^3.  The conventional
hexagonal reference parameters are `a=5.0355 A` and `c=13.7471 A`.

| Folder | V/Vref | V primitive (A^3) | a hex (A) | c hex (A) |
|---|---:|---:|---:|---:|
| V092 | 0.92 | 92.574952 | 4.897471 | 13.370276 |
| V094 | 0.94 | 94.587451 | 4.932706 | 13.466468 |
| V096 | 0.96 | 96.599950 | 4.967444 | 13.561306 |
| V098 | 0.98 | 98.612449 | 5.001704 | 13.654835 |
| V100 | 1.00 | 100.624948 | 5.035500 | 13.747100 |
| V102 | 1.02 | 102.637447 | 5.068849 | 13.838143 |
| V104 | 1.04 | 104.649946 | 5.101764 | 13.928004 |
| V106 | 1.06 | 106.662445 | 5.134260 | 14.016719 |
| V108 | 1.08 | 108.674944 | 5.166350 | 14.104326 |

The QE `celldm(1)` values use the corresponding isotropic linear scaling
factor `(V/Vref)^(1/3)`.  `celldm(4)` is unchanged.

## Run the nine calculations

Enter this directory on Jakar and submit the Slurm array:

```bash
cd Fe2O3/eos_birch_murnaghan
sbatch run_eos_array.sbatch
```

This submits array tasks 0 through 8, with at most three running
simultaneously.  Each task enters one volume folder and writes
`alpha_Fe2O3.eos.out` there.  The inputs find the pseudopotentials in the
parent `Fe2O3/pseudo/` folder through `../../pseudo`.

Check completion with:

```bash
grep -l 'JOB DONE' V*/alpha_Fe2O3.eos.out
```

There must be nine returned files.  Also check the final magnetic moments and
forces:

```bash
for d in V0*; do
  echo "$d"
  grep 'Total force' "$d/alpha_Fe2O3.eos.out" | tail -1
  grep -E 'total magnetization|absolute magnetization' "$d/alpha_Fe2O3.eos.out" | tail -2
done
```

Every point must remain in the same AFM electronic state.  A collapsed or
different magnetic solution produces an invalid discontinuity in the EOS.

## Fit the data

The fitter can also produce a preliminary result after at least five jobs
finish:

```bash
python3 fit_birch_murnaghan.py
```

Incomplete points are skipped.  For `V100`, the script first checks
`V100/alpha_Fe2O3.eos.out`; if that calculation is incomplete, it
automatically tries the completed parent phonon relaxation at
`../alpha_Fe2O3.relax.out`.  The fallback is accepted only when its actual
volume matches the `V100` target within 0.5 percent.

The script extracts the final QE energy and actual cell volume and fits the
third-order Birch-Murnaghan equation.  It creates:

- `eos_data.csv`
- `eos_results.txt`
- `alpha_Fe2O3_birch_murnaghan.png`
- `alpha_Fe2O3_birch_murnaghan.pdf`

SciPy is used when available.  A slower NumPy-only fitting method is included
as a fallback.  Plotting requires matplotlib but the numerical fit does not.

Inspect the plot rather than accepting the fit blindly.  The minimum should
be bracketed by the calculated volumes, the curve should be smooth, and the
fit residual should be small.  If the fitted minimum lies near either end,
add more points on that side before reporting `V0`, `B0`, or `B0'`.

A five-point partial fit is explicitly labeled `PRELIMINARY`.  It is useful
for checking the trend, but it should not replace the final nine-point fit,
especially if all completed points lie on only one side of the minimum.

For final publication values, also converge the k grid and repeat the EOS for
the selected Hubbard U.  The fitted volume and bulk modulus can depend
noticeably on U.
