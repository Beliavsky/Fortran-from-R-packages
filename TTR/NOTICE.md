# License, attribution and provenance

This is a derivative computational translation of **TTR 0.24.4**, “Technical Trading Rules,” into modern free-form Fortran. It is distributed under **GNU GPL version 2 or, at your option, any later version**. The full GPL version 2 is in `LICENSE`. No warranty is provided; in particular there is no warranty of merchantability or fitness for a particular purpose.

## Upstream

- Author and maintainer: **Joshua M. Ulrich**.
- Contributing author listed in DESCRIPTION: **Ethan B. Smith**.
- Source: https://cran.r-project.org/src/contrib/Archive/TTR/TTR_0.24.4.tar.gz
- Project: https://github.com/joshuaulrich/TTR
- CRAN publication date: 2023-11-28.
- Official source archive SHA-256: `89732b9c359bae2f41cd23db649f0897c10fab0702d780c4c25a997322710284`.

Original R and C computational/reference source files, original DESCRIPTION, NAMESPACE and THANKS are retained unmodified in `upstream/` for algorithm auditing and preservation of original notices. They are not compiled or required at runtime. `upstream/SHA256SUMS` records their individual hashes. These are the translated package's original sources, not vendored dependency implementations.

Applicable source notices include:

- Copyright (C) 2007–2013, 2007–2017, 2007–2018, 2017 and 2020 **Joshua M. Ulrich**, as individually stated in upstream source headers.
- Copyright (C) 2012–2017 **Charlie Friedemann and Joshua M. Ulrich**, for percent-rank code.
- Copyright (C) 2016 **Peter Carl and Joshua M. Ulrich**, for SNR.
- Price Bands implementation by **Brian G. Peterson**, inspired by univariate Bollinger Bands and **Ram Ben-David**.
- DVI indicator attribution to **David Varadi**.
- Upstream THANKS credits **Jeff Ryan**, including his example EMA Fortran implementation, and **Ion Georgiadis** for suggestions and testing.

All other source-specific algorithm references and attributions remain intact in the retained upstream files. New maintained Fortran, examples and test scaffolding in this derivative are provided under GPL-2.0-or-later. Translation and verification changes were made in 2026; this is not an upstream TTR release.

## Dependency reuse

The repository was checked before implementation at revision `cbd437d6d92f5481706fd8ea5a7829aaff7b6e24`. Existing `rfortran-core` is used through `path = "../rfortran-core"`. Its single `dp` kind, descriptive statistics, median and sorting routines are reused where compatible. It remains separately licensed under MIT with its original authors/notices in that sibling package; no dependency source is copied here.

Existing shared rolling helpers were inspected. TTR's leading-missing validation, cumulative startup, exact C rolling-sum update order, unusual Wilder seed, EMA/ZLEMA initialization, asymmetric percent-rank tie comparison and stateful rules require package-specific code around compatible shared statistical operations. In particular rolling update order affects downstream rank ties, so a merely algebraically equivalent shared sum is insufficient. No BLAS, LAPACK, ARPACK, r.f90, r_mod.f90 or translated dependency sources are embedded.

## Oracle provenance

Reference generators call installed, unmodified TTR 0.24.4 using R 4.5.0 in the validation environment. Oracle dependencies were zoo 1.9-1, xts 0.14.3, curl 8.0.0 and lattice 0.23-1. The checked-in Fortran expected values and text data are deterministic generated test data. No R runtime is required for FPM tests. Two repaired Yang–Zhang override cases use the published formula with genuine upstream primitives because the corresponding upstream wrapper branch errors; they are identified in the generator.
