# Porting notes

## Numerical design

The original R package dispatches most calculations to the header-only C++
`vinecopulib` library. This translation implements the continuous parametric
algorithms directly in Fortran.

Archimedean families share generator, inverse-generator, and derivative code.
The copula density is evaluated from the standard generator derivative
identity. Gaussian and Student-t h-functions are analytic. Their bivariate CDFs
are obtained by integrating the conditional CDF. Tawn uses its Pickands
function and derivatives.

All probability arguments are clipped away from exactly zero and one to avoid
infinite normal or Student quantiles. Inverse h-functions use monotone
bisection, which is slower than specialized closed forms but stable across all
families and rotations.

## Vine representation

A C-vine pair `(tree, edge)` is stored as `pair(tree, edge)`, where the tree
root is `tree` and `edge > tree`. A D-vine pair for variables at ordered
positions `i` and `j`, conditional on the variables between them, is stored as
`pair(i, j)`.

The Rosenblatt recursions use conditional CDF tables. Inverse Rosenblatt
transforms solve each sequential conditional equation by bisection. This is
simple and reliable, but slower than the optimized triangular-array C++
implementation for high dimensions.

## Estimation

Bivariate estimation uses bounded Nelder-Mead. Starting values are based on
sample Kendall's tau for common families and conservative defaults for
multi-parameter families. Family selection can use AIC, BIC, or log
likelihood. C-vine and D-vine estimation proceeds sequentially from lower to
higher trees.

## Unsupported upstream features

The upstream package's strongest performance advantages come from optimized
C++ code, arbitrary R-vine structures, nonparametric copulas, discrete-variable
support, multithreading, and structure selection. Those are not represented as
completed features in this native port. The retained upstream snapshot is
included so future work can be checked against the exact attached source.
