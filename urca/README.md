# urca-fortran

Modern Fortran 2018 translation of the computational core of the R package
`urca` 1.3-4 (Pfaff, Zivot, and Stigler).

The port is intended for direct use from Fortran and as a numerical library
that can later be wrapped from Python, R, Julia, MATLAB/Octave, or C.
It uses FPM for project layout and BLAS/LAPACK for dense linear algebra.

## Implemented statistical functionality

### Univariate unit-root and stationarity tests

- Augmented Dickey-Fuller (`ur.df`)
  - none, drift, and trend specifications
  - fixed, AIC, and BIC lag selection
  - tau and upstream phi statistics/critical values
- Elliott-Rothenberg-Stock (`ur.ers`)
  - DF-GLS and P-test
  - constant and trend specifications
- KPSS (`ur.kpss`)
  - level and trend stationarity
  - short/long or explicit HAC lag choices
- Phillips-Perron (`ur.pp`)
  - Z-alpha and Z-tau
  - constant and trend specifications
  - auxiliary Z-tau-mu/Z-tau-beta statistics
- Schmidt-Phillips (`ur.sp`)
  - tau and rho statistics
  - polynomial degrees 1 through 4
- Zivot-Andrews (`ur.za`)
  - intercept, trend, and both-break models
  - break-date search and lagged differences

### Cointegration and VECM calculations

- Phillips-Ouliaris (`ca.po`)
  - Pu and Pz statistics
  - none/constant/trend deterministic terms
- Johansen (`ca.jo`)
  - trace and maximum-eigenvalue tests
  - none/constant/trend deterministic terms
  - long-run and transitory specifications
  - seasonal dummies and user-supplied deterministic regressors
  - eigenvectors, loading matrices, residual covariance matrices, and VECM
    coefficient objects
- `cajools`, `alphaols`, and `cajorls`
- Level-shift Johansen procedure (`cajolst`)
  - endogenous break-date search
  - shifted deterministic terms
  - dedicated critical-value tables

### Johansen restrictions

Direct numerical counterparts of:

- `blrtest`
- `alrtest`
- `ablrtest`
- `bh5lrtest`
- `bh6lrtest`
- `lttest`

The partly-known-beta procedures retain the iterative algorithms used by the
upstream package rather than replacing them with generic Wald tests.

### MacKinnon response surfaces

- `punitroot`
- `qunitroot`
- `unitrootTable`

The response-surface coefficient data embedded in `urca` are translated to
Fortran parameter arrays. The numerical GLS interpolation/extrapolation
algorithm follows the bundled `UnitRootMacKinnon.f` implementation.

## Not translated

The following are intentionally outside the numerical-library scope:

- R S4 class definitions and method dispatch
- R formula/model-frame construction
- `show`, `summary`, and table formatting
- plotting (`plot`, `plotres`)
- R-specific time-series indexing and `na.omit` behavior
- bundled `.rda` datasets and book examples as runtime objects

The original package tree is retained under `upstream/` for provenance and
comparison.

## nlme

The original `urca` NAMESPACE imports `nlme`, and the user supplied a modern
Fortran translation of `nlme`. A source audit found no call from the exported
`urca` computational routines to an `nlme` numerical routine. Consequently,
`nlme-fortran` is retained under `reference/nlme-fortran/` but is deliberately
not linked as a dependency. This keeps the numerical port smaller and avoids
adding a dependency that the translated algorithms do not use.

## Building with FPM

A BLAS/LAPACK installation is required.

```text
fpm build
fpm test
fpm run --example urca_example
```

The manifest links `lapack` followed by `blas`.

The source uses standard free-form line lengths and does not require GNU
Fortran's `-ffree-line-length-none` extension.

## Direct GNU Fortran validation

The release was validated with GNU Fortran 14.2.0 using Fortran 2018,
runtime checking, and implicit interfaces promoted to errors. For example:

```text
gfortran -std=f2018 -fcheck=all -Werror=implicit-interface ... -llapack -lblas
```

See `VALIDATION.md` for the regression references.

## License and provenance

The upstream `urca` package is licensed under GPL version 2 or, at the user's
option, any later version. The translated project therefore uses
`GPL-2.0-or-later`.

The MacKinnon response-surface material has additional provenance documented
by the upstream package. In particular, `upstream/inst/Licenses/` contains the
upstream license text and the 2008 permission correspondence concerning use
of the MacKinnon code in `urca` under the GPL. These files are retained
unchanged.

See `LICENSE`, `NOTICE.md`, and `PORTING_NOTES.md`.
