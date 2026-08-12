# Upstream provenance

Source archive: `lpSolve-main.zip`

Upstream DESCRIPTION:

- Package: `lpSolve`
- Version: `5.6.23.9000`
- Title: Interface to `lp_solve` v5.5 to Solve Linear/Integer Programs
- Authors: Michel Berkelaar; current R package maintenance by Gabor Csardi
- License declaration: `LGPL-2`

The complete supplied upstream tree is retained under `original/lpSolve-main/`.
The original tree contains the lp_solve 5.5 C implementation and third-party
support components such as LUSOL/COLAMD; none of those C files are compiled by
the FPM target.
