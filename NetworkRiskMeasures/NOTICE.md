# Upstream provenance, license and dependency reuse

This is a computational translation of **NetworkRiskMeasures 0.1.7** by
**Carlos Cinelli** (author, maintainer) and **Thiago Cristiano Silva** (author).
The package DESCRIPTION declares **GPL-3**. The complete GNU General Public
License version 3 is included as LICENSE; the translation is GPL-3.0-only.
No upstream copyright year is invented; original authorship and metadata are
retained verbatim in upstream/DESCRIPTION and the source comments.

- CRAN: https://cran.r-project.org/package=NetworkRiskMeasures
- Official source: https://cran.r-project.org/src/contrib/NetworkRiskMeasures_0.1.7.tar.gz
- Source archive SHA-256: 245c899f97678e3716af0b96649c03412b5a6f1ae474b61bb106eeaffc3a578d
- Upstream repository: https://github.com/carloscinelli/NetworkRiskMeasures
- CRAN publication: 2025-05-04 19:20:02 UTC
- Target: https://github.com/Beliavsky/Fortran-from-R-packages
- Dependency inspection at repository commit: cbd437d6d92f5481706fd8ea5a7829aaff7b6e24

upstream/DESCRIPTION, upstream/NAMESPACE and upstream/R/*.R preserve the exact
original metadata, computational source, documentation, references and notices.
upstream/SHA256SUMS records their hashes. Source files are reference material,
not compiled. R-specific containers, plotting, printing, interactive behavior,
data-frame construction, sparse storage and datasets are omitted from the
Fortran implementation. The numerical body of registered summary.contagion is
also translated as a helper within the contagion mapping.

## Research attribution preserved in upstream source

- Estrada and Hatano (2008), communicability in complex networks,
  Physical Review E 77:036111.
- Silva, Souza and Tabak (2015), Central Bank of Brazil Working Papers 391/392,
  network structure, monitoring vulnerability and impact diffusion.
- Upper and Worm (2004), estimating bilateral interbank exposures,
  European Economic Review 48, 827–849.
- Anand, Craig and von Peter (2015), filling in the blanks: network structure
  and interbank contagion, Quantitative Finance 15(4), 625–636.
- Bardoscia, Battiston, Caccioli and Caldarelli (2015), DebtRank microscopic
  foundation, PLoS ONE 10(6):e0130406, doi:10.1371/journal.pone.0130406.
- Battiston, Puliga, Kaushik, Tasca and Caldarelli (2012), original DebtRank,
  Scientific Reports 2:541.

## Reused implementations, never vendored

Existing top-level translations were inspected before implementation:

- rfortran-core supplies r_kinds:dp. Its manifest identifies Vivek Rao and
  contributors and MIT licensing.
- expm (FPM name expm-fortran) supplies expm_module:expm. Its unchanged
  GPL-3.0-or-later manifest uses sibling rfortran-linalg. It implements balanced
  Higham scaling/squaring and Padé algorithms rather than a new local copy.
- rfortran-linalg supplies the shared numerical infrastructure transitively.
  Its unchanged manifest pins https://github.com/perazz/fortran-lapack.git at
  9982df3b660918761b993eaf0ff2e507547e25b2. **Federico Perini**, the original
  fortran-lapack author, and its LAPACK/BLAS contributors retain attribution
  and licenses in that dependency. This package makes no direct LAPACK calls
  and introduces no direct LAPACK dependency. The requested Beliavsky fork
  (https://github.com/Beliavsky/fortran-lapack) would be the choice for new
  direct LAPACK access; the compatible existing shared package is reused
  unchanged here rather than rewriting its dependency manifest.

No dependency source, r.f90, r_mod.f90, BLAS, LAPACK or ARPACK source is copied
into this package or its release ZIP. No shared package has been edited.
The result builds using explicit module interfaces and pure Fortran dependency
sources, without external numerical-library linking. Communicability wrappers
are impure solely because the unchanged expm API is impure.

The Fortran port changes types, validation, diagnostics, storage and RNG as
detailed in API.md. It preserves directed orientation, weighted definitions,
raw contagion stress, minimum-density objective quirks and convergence tests.
See README.md and fpm.toml for the identical 12/12 coverage accounting.
