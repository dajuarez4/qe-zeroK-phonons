# TDEP 6 Angstrom cutoff test

This directory is a separate 6.0 Angstrom second-order-cutoff test using the
same 413 complete MD configurations as the main 5.5 Angstrom TDEP fit. The
unfinished 414th coordinate block was excluded.

## Result

- Requested second-order cutoff: 6.0 Angstrom
- TDEP-reported maximum cutoff for this supercell: 5.67126 Angstrom
- Complete configurations: 413
- Represented trajectory time: 399.6188 fs
- Harmonic force-fit residual RMS: 0.330090 eV/Angstrom
- Force residual R-squared: 0.901834
- Anharmonicity measure: 0.313315
- Frequency range: -1.439419 to 25.897203 THz

The fitted force constants and dispersion are byte-for-byte identical to the
5.5 Angstrom results. Thus, this supercell contains no additional independent
interaction shell between the 5.5 Angstrom model and its 5.67126 Angstrom
geometric reach; requesting 6.0 Angstrom does not enlarge the fitted model.
