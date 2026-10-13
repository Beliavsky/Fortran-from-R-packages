# License and provenance

This package translates the computational routines of R package **irr 0.85**
(release date 2026-06-17), by Matthias Gamer, Jim Lemon, Ian Fellows and
Puspendra Singh. Upstream license: **GNU General Public License version 2 or any
later version**. This translation and the retained R sources use that same
license. See the complete GPL version 2 text in COPYING. There is no warranty.

Upstream sources:

- CRAN: https://cran.r-project.org/package=irr
- Source release: https://cran.r-project.org/src/contrib/irr_0.85.tar.gz
- CRAN mirror: https://github.com/cran/irr
- Inspected mirror commit: `98e4535cece091e099e7ad7ad96fb132a53fc2d8`
- Release archive SHA-256:
  `f32e5a703806c5a48ac95997ee1f4abfe0df47f0c7cecc19a5c35be9908ed42a`

The original DESCRIPTION, NAMESPACE, CHANGES and full R source tree are retained
under `upstream/`, without rewriting source comments or author attribution.
`upstream/SHA256SUMS` identifies the preserved files. The two printing functions
are retained for provenance but do not contribute to translation coverage.

The original comments cite Finn (1970), Bartko (1966), Shrout and Fleiss (1979),
McGraw and Wong (1996), Asendorpf and Wallbott, Janson and Olsson (2001), and
Eliasziw et al. (1994) as the methodological sources. Their original references
and explanatory comments remain in the retained R sources. R documentation
used for fixtures and arguments is retained under `upstream/man/`.

## Changes in this translation

The 2026 Fortran translation replaces R object handling with explicit typed
arrays, derived result types and status codes; adds shape/finite/resource
checks; computes upper tails directly; and makes pure/elemental procedures
where the dependency call chain allows. The LP path remains impure. Behavioral
differences, preserved upstream quirks and optional mathematical corrections
are described in README.md. No authorship of the original algorithms is claimed.

## Dependencies, not bundled

- `rfortran-core`: shared numerical/statistical helpers, MIT, by Vivek Rao and
  contributors. `dp`, ordering, correlations and distribution routines are reused.
- `rfortran-linalg`: shared high-level linear algebra, MIT, by Vivek Rao and
  contributors. The translation uses its explicit pure `solve_system` API.
- `lpSolve`: sibling `lpsolve-fortran` translation, LGPL-2.0-only, using its
  continuous LP solver for N2 variance maximization. No solver source is copied.

These siblings were tested from Beliavsky/Fortran-from-R-packages revision
`cbd437d6d92f5481706fd8ea5a7829aaff7b6e24`. Their license notices remain with those
packages. The tested shared linalg manifest pins Federico Perini's
`perazz/fortran-lapack` at `9982df3b660918761b993eaf0ff2e507547e25b2`.
It supplies Fortran BLAS/LAPACK without separate system libraries. The translation
needs no direct LAPACK access; therefore it does not override that tested shared
manifest. The FPM-compatible Beliavsky fork
https://github.com/Beliavsky/fortran-lapack is the requested route for direct
LAPACK or other pure chains, but is **not used directly by this package**.
No dependency code, BLAS/LAPACK source, or copied R compatibility module is
included in the deliverable.
