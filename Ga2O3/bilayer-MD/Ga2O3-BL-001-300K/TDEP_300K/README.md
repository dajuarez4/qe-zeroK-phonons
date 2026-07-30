# Updated diagnostic TDEP phonons from the partial QE trajectory

TDEP completed a second-order force-constant fit and a phonon dispersion
calculation along the reciprocal-fractional path
`Gamma-X-S-Y-Gamma`.

This result is a pipeline test, **not a physically converged 300 K phonon
spectrum**.

## Data actually available

- Complete MD configurations: 20
- Simulated time represented: 19.352 fs
- Mean sampled temperature: 417.740 K
- Temperature of latest included frame: 616.859 K
- Incomplete configurations excluded: 1
- Simulation cell: the 10-atom primitive bilayer cell

The QE output currently ends after printing the position and temperature
for step 21 but before printing its matching force. That frame was
excluded. Its reported temperature is 646.103 K.

## Fit diagnostics

- Force equations: 600
- First- plus second-order fit parameters: 135
- Overdetermination ratio: 4.4
- Cross-validation energy R2: 0.425
- Total residual-force RMSE: approximately 0.69 eV/Angstrom
- Dispersion range: approximately -21.9 to 27.7 THz

The many large negative frequencies are not evidence of a physical
instability. The trajectory is heating rapidly instead of equilibrating
near 300 K, and the primitive simulation cell cannot resolve spatial force
constants. TDEP also reports a maximum unaliased cutoff of only 1.571
Angstrom, smaller than the nearest interaction distance, while this
diagnostic fit requested 2.5 Angstrom.

Ga2O3 is polar, but this run contains no Born effective charges or
dielectric tensor, so non-analytic LO-TO corrections are absent.

## Files

- `outfile.forceconstant`: fitted diagnostic second-order force constants
- `outfile.dispersion_relations`: 30 branches in THz
- `Ga2O3-BL-001-TDEP-diagnostic.png`: diagnostic dispersion plot
- `extract_forceconstants.log`: fit settings and diagnostics
- `phonon_dispersion_relations.log`: dispersion calculation log

For a production TDEP result, first relax the bilayer, run MD in an
in-plane supercell large enough for the desired cutoff (approximately
`4 x 2 x 1` for a 5 Angstrom target), equilibrate at 300 K, and collect
hundreds of decorrelated configurations.
