# Notices, licensing and upstream provenance

## Upstream

- R package: Kendall, version 2.2.2
- Author and maintainer: A.I. McLeod, as recorded in upstream/DESCRIPTION
- License: GPL (>= 2), expressed here as GPL-2.0-or-later
- CRAN mirror repository: https://github.com/cran/Kendall
- Exact translated revision: 25bac7757bfdfbae4d47f229b68ccb2db64bc3b4
- Immutable source tree:
  https://github.com/cran/Kendall/tree/25bac7757bfdfbae4d47f229b68ccb2db64bc3b4
- Target repository was inspected at cbd437d6d92f5481706fd8ea5a7829aaff7b6e24:
  https://github.com/Beliavsky/Fortran-from-R-packages

The upstream DESCRIPTION, NAMESPACE, NEWS, all R sources, three computational
manual pages and native source are preserved byte-for-byte under upstream/.
For build isolation only, original `src/Kendall.f` is stored as
`upstream/src/Kendall.f.txt`; its SHA-256 is
`21294db32acc7f30fe20fba1c57f1e594d87f6c32d0cb84f449c9445895e0287`.
Other paths retain their upstream basenames. SHA256SUMS lists all preserved
files. Dataset binaries and presentation-only manuals are not included.
The upstream DESCRIPTION's date and publication metadata are preserved as supplied.

The maintained source `src/kendall_module.f90` translates the computational
formulas in `R/Kendall.R`, `R/MannKendall.R`, `R/SeasonalMannKendall.R` and
`src/Kendall.f`. The modern port uses explicit interfaces, pure procedures,
binary64 arithmetic, typed results, explicit statuses, shared helpers, and
an independent optional no-ties exact distribution recurrence. Changes were
prepared in 2026 for the target repository. They are distributed under
GPL-2.0-or-later, like upstream. No additional upstream copyright years or
copyright holders have been invented.

The upstream native source attributes Algorithm AS 66 (normal integral,
Hill 1973) and AS 71 (Kendall upper-tail probability, Best and Gipps 1974).
Their original comments remain preserved. Neither their native implementation
nor the unreachable AS 71 branch is used by the maintained library:
normal tails use the shared dependency and the optional exact method uses a
nonnegative inversion-probability recurrence. The Kendall tie formulas follow
M.G. Kendall, Rank Correlation Methods; the seasonal aggregation follows the
upstream R wrapper. Original bibliographical references remain in upstream/man/.

## Shared dependency reuse

`rfortran-core` is an existing sibling package in the target repository:
`rfortran-core = { path = "../rfortran-core" }`.
It is MIT-licensed and credits Vivek Rao and contributors in its manifest.
The dependency remains separate and retains its own LICENSE. This package imports:

- `r_kinds: dp`, the single shared real kind
- `r_ordering: r_sort_values_in_place`, stable merge sorting for tie groups
- `r_distributions: r_pnorm`, a direct upper-tail normal probability

The repository's shared package was checked before implementation. The existing
ppcor translation contains a private tau-b helper but does not expose a compatible
pure, scalar Kendall result API with scores, variances and p-values. It is not
copied. No BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, dependency source, or legacy
compatibility package is embedded here. No system numerical libraries are linked.

## License text

LICENSE contains the unmodified GNU General Public License version 2 text.
The upstream grant explicitly permits version 2 or any later version; including
the version 2 text does not narrow that grant.
