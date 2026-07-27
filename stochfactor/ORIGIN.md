# Origin and licensing

This project was translated from the two attached source archives:

- `stochvol-master.zip`
  - Package: `stochvol`
  - Version: 3.2.9
  - Declared license: `GPL (>= 2)`
  - Archive SHA-256: `e27fb2d7d19e5115573c8fcc96940d12f89f9e28e450c34d5b66dab96bf6275a`
  - Source archive commit marker: `b05b04bd3b76a08e546ec15f334832bb74b5d8f1`
  - Original authors listed in DESCRIPTION: Darjus Hosszejni and Gregor Kastner

- `factorstochvol-master.zip`
  - Package: `factorstochvol`
  - Version: 1.1.2
  - Declared license: `GPL (>= 2)`
  - Archive SHA-256: `332fa8ad170ed03256d8ed997e6193f3f221823286f9df3afb353b1a2d73a640`
  - Source archive commit marker: `648b212e8e8ec942ff183589c2cd765007fe4472`
  - Original author/contributors listed in DESCRIPTION: Gregor Kastner,
    Darjus Hosszejni, and Luis Gruber

`factorstochvol` imports and links to `stochvol`. The Fortran translation preserves
that dependency structurally by placing the reusable univariate SV machinery in
`sv_*` modules and calling it from `fsv_core`.

Both packages allow redistribution under GNU GPL version 2 or any later version.
The combined project uses `SPDX-License-Identifier: GPL-2.0-or-later` in every
Fortran source and test file. `LICENSE` contains the complete GNU GPL version 2
text.

The Fortran code is a clean-language translation and numerical reimplementation.
It is not endorsed by the original package authors.
