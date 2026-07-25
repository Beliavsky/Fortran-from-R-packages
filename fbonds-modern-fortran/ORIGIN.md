# Origin and Provenance

This project is a modern Fortran translation of the computational content of:

```text
Package: fBonds
Version: 3042.78
Title: Rmetrics - Pricing and Evaluating Bonds
Date: 2017-11-12
Authors: Diethelm Wuertz, Tobias Setz
License: GPL (>= 2)
```

The translated source archive was `fBonds-master.zip` supplied by the user.
The package contains one numerical R source file, `R/TermStructure.R`, exporting
`NelsonSiegel` and `Svensson`.

## Licensing

The original metadata specifies `GPL (>= 2)`. The translation therefore uses:

```text
SPDX-License-Identifier: GPL-2.0-or-later
```

`LICENSE` contains the GNU General Public License version 2 text. Every Fortran
source, application, example, and test file carries the GPL-2.0-or-later SPDX
identifier and notice.

## Documented source details

### Nelson-Siegel starting grid

The R grid constructs normal equations from the standard Nelson-Siegel yield
loadings, but its temporary `yfit` line uses an exponential basis. The final
curve and optimizer use the standard yield loadings. The Fortran implementation
uses the standard loadings consistently in the grid and final fit.

### Svensson objective

The R objective function evaluates an SSE expression and then an L1 expression.
Because R returns the final expression, `nlminb` effectively minimizes L1. The
Fortran default preserves L1 and also exposes an explicit SSE option.

### Optimizer

R's `nlminb` is replaced by bounded Nelder-Mead. Positive decay parameters are
optimized on the log scale. This is a tested numerical analogue; exact iteration
or endpoint equivalence is not claimed.
