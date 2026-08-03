# Porting notes

## Interface model

R vectors, matrices, lists, closures, missing values, and dynamically dispatched
functions were replaced by explicit Fortran arrays, procedure interfaces,
derived result types, optional arguments, and integer status codes.  R names
containing periods generally use underscores in Fortran.  The aggregate
`pracma` module re-exports all translated modules.

## Important numerical adaptations

- **Linear algebra:** dense, self-contained algorithms are used.  Symmetric
  eigenproblems are the primary supported path for `eig`; singular-value and
  pseudoinverse operations are derived from symmetric eigensystems.
- **FFT:** the portable implementation is a direct discrete Fourier transform.
  It preserves results but not the asymptotic speed of a vendor FFT library.
- **Optimization:** portable line-search, Nelder-Mead, pattern-search,
  least-squares, penalty, linear-programming, NNLS, and Goldfarb-Idnani QP
  routines replace R package callbacks and external optimizers.
- **ODE aliases:** RK4, ABM3, RKF45, Newmark, and Crank-Nicolson are direct
  implementations.  Several higher-order R solver names share the adaptive
  RKF45 backend rather than reproducing every original Butcher tableau.
- **Interpolation:** interpolation objects are typed derived structures.
  Closure-returning R helpers such as `pchipfun` are represented through
  coefficients plus `ppval`, not dynamically generated functions.
- **Randomness:** deterministic portable generators are supplied.  Streams do
  not match R's RNG bit for bit.
- **Magic squares:** standard odd-order and doubly-even constructions are
  implemented; difficult even-order compatibility paths use a documented
  fallback.
- **Boundary-value problems:** `bvp` supplies a finite-difference scalar
  second-order boundary-value solver.  R's general shooting/closure framework
  is not reproduced in full.
- **Kriging and spatial helpers:** portable dense formulations are intended for
  moderate problem sizes.

## Missing values

Where meaningful, IEEE NaNs are recognized as missing/nonfinite values.  The
handling policy is procedure-specific and is documented in `API.md`; unlike R,
there is no universal `NA` object or automatic recycling.

## Deliberately omitted R-only functionality

The compiled library omits:

- plotting and graphics helpers (`figure`, `errorbar`, `ezplot`, `ezcontour`,
  `ezsurf`, `ezmesh`, `polar`, `quiver`, `plotyy`, `semilog*`, `loglog`, and
  `andrewsplot`);
- workspace, filesystem, timing, and console helpers (`clear`, `who`, `whos`,
  `what`, `cd`, `pwd`, `tic`, `toc`, `beep`, `disp`, `fprintf`, `ver`, and
  `matlab`);
- R regular-expression and string-environment helpers (`regexp`, `regexpi`,
  `regexprep`, `refindall`, `strfind`, `strcat`, `strcmp`, `num2str`,
  `str2num`, and related routines);
- dynamic closure factories and formula/data-frame behavior.

## Error handling

Invalid dimensions, impossible parameter combinations, singular systems, and
nonconvergence return explicit status codes and, where appropriate, zero-size
arrays or IEEE NaNs.  The library avoids terminating the caller for ordinary
numerical failures.
