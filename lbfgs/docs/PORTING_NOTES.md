# Porting notes

## Source lineage

The R package wraps Naoaki Okazaki's libLBFGS implementation. The modern
Fortran solver was translated from the supplied package's `src/lbfgs.cc` and
`src/lbfgs.h`, including the More-Thuente and backtracking line searches and
the OWL-QN pseudo-gradient and orthant projection operations.

No external Fortran optimization library is required. Original files are
retained in `upstream/lbfgs-master`.

## Intentional Fortran adaptations

- Arrays are ordinary allocatable Fortran arrays; SSE-specific aligned
  allocation and padded dimensions are unnecessary.
- OWL-QN index ranges use one-based inclusive indexing.
- Callback state uses optional unlimited-polymorphic `user_data` rather than
  an R environment or `void *`.
- The objective and gradient can be evaluated together, avoiding duplicate
  work. A separate-callback overload remains available for direct R API
  correspondence.
- Non-finite initial evaluations terminate with `lbfgserr_logicerror`.
  Non-finite trial evaluations in backtracking are rejected by shortening the
  step.
- Degenerate interpolation denominators are safeguarded with midpoint or
  boundary fallbacks.
- A nonpositive curvature pair is skipped rather than used to divide by zero.

## Status-code discrepancy in the R wrapper

The supplied `src/lbfgs.h` uses libLBFGS status values beginning at `-1024`.
The R file's message-selection chain checks values beginning at `-1204`.
The Fortran port follows the actual numerical library constants in the header,
not the inconsistent R message checks.

## Defaults

The Fortran defaults follow the exported R function, including
`max_linesearch = 20` and `delta = 0`, rather than libLBFGS's lower-level C
structure defaults where those differ.
