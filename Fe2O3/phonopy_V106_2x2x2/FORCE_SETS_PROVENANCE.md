# V106 FORCE_SETS provenance

`disp-003/alpha_Fe2O3.fd.scf.out` reached `JOB DONE` but did not converge
after 300 electronic iterations (final estimated SCF error: 0.24718016 Ry).
Consequently, QE did not print a force block and Phonopy could not parse it.

The current `FORCE_SETS` is a harmonic reconstruction made on 2026-08-19.
Displacements 003 and 004 are an exact positive/negative pair. The missing
force was reconstructed as

    F003 = 2 F0 - F004,

where `F0` is the average residual-force baseline from the two complete
central pairs, `(F001 + F002 + F005 + F006) / 4`, after Phonopy's drift-force
subtraction. The baseline RMS was about 8.9e-6 in the QE force units when
expressed as the correction to the simple `F003 = -F004` estimate; the RMS
force in displacement 004 was 3.77e-4 in the same units.

The resulting bands, DOS, and thermal properties are suitable as a harmonic
estimate. For a fully ab-initio V106 dataset, rerun displacement 003 to true
SCF convergence and regenerate all derived files with `collect_forces.sh`.
