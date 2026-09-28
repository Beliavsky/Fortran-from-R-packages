# API coverage

## Strict R-function translation coverage

The coverage denominator follows the repository convention: exported functions
that merely dispatch to another R package are excluded. In rARPACK 0.11-0 all
three exports are one-line calls into RSpectra, so there are **0 exported
computational R functions in the strict coverage basis**. The manifest records
0 translated of 0 total and uses `0.0` as the numeric fraction convention for
this empty denominator.

## Fortran compatibility API

Although the wrapper-only R exports are excluded from the strict translation
fraction, the Fortran package provides useful numerical counterparts:

| R-facing operation | Fortran procedure | Coverage | Material differences |
| --- | --- | --- | --- |
| `eigs` | `rarpack:eigs` | substantial | Dense real matrices only; returns complex eigenvalues/vectors; supports LM/SM/LR/SR/LI/SI and real shift-and-invert. R `Matrix` sparse classes and R function callbacks are omitted. |
| `eigs_sym` | `rarpack:eigs_sym` | substantial | Dense real matrices; honors lower/upper triangle; supports LM/SM/LA/SA/BE and real shift-and-invert. R sparse classes and callbacks are omitted. |
| `svds` | `rarpack:svds` | substantial | Dense real matrices; leading singular triplets and independent `nu`/`nv` counts are supported. The implementation uses ARPACK on the smaller normal-equation operator rather than RSpectra internals. |

The R `opts` list is represented by explicit optional Fortran arguments `ncv`,
`tol`, `maxitr`, and for `eigs`/`eigs_sym`, `retvec`. Result derived types carry
the numerical outputs plus `nconv`, `niter`, `nops`, and `info` status fields.

## Status codes

ARPACK status codes are returned unchanged where applicable. Package-level
argument validation uses negative codes below -1000; callers should treat any
nonzero `info` as a non-success status.
