# Notices and attribution

This project is a computational translation of **matrixdist 1.1.9**, originally authored by Martin Bladt and Jorge Yslas, with Alaric Mueller listed as contributor.  Upstream declares `License: GPL-3` and is preserved here as GPL-3.0-only for translated matrixdist code.

The supplied `r_mod.f90` is separate MIT-licensed code and retains its MIT notice.

The upstream DESCRIPTION identifies the numerical methodology with, among others, the following work: Asmussen, Nerman & Olsson (1996) on EM fitting of phase-type distributions; Olsson (1996) on censored phase-type estimation; Albrecher & Bladt (2019); Albrecher, Bladt & Yslas (2022); Bladt & Yslas (2022/2023); Bladt (2022/2023); and Albrecher, Bladt & Mueller (2023).  The original DESCRIPTION is retained under `upstream/` with the DOI information exactly as supplied by the package.

The Fortran port preserves the package's Van Loan matrix-integral approach and phase-type EM sufficient-statistic formulas.  Original R and C++ computational sources are retained under `upstream/R` and `upstream/src` for auditability and attribution.
