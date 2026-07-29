# Notices and provenance

This project is a modern Fortran translation of ufRisk 1.0.7.

Original ufRisk authors listed by the package metadata:

- Yuanhua Feng
- Xuehai Zhang
- Christian Peitz
- Dominik Schulz
- Shujie Li
- Sebastian Letmathe

The original ufRisk source and metadata are retained under `original/ufRisk/`.

The user-supplied rugarch Fortran inputs are retained unmodified under `original/rugarch-input/`; the integrated copies in `src/` contain small build-safety and integration changes documented in `PORTING.md`.

The library also incorporates source from the earlier modern Fortran translations of `fracdiff` and `smoots`. Their source headers retain their original SPDX identifiers. The combined work is distributed under GPL-3.0-only, which is compatible with the incorporated GPL-2.0-or-later and GPL-3.0-only components.

No endorsement by the original R-package authors is implied.
