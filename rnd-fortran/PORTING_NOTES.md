# Porting notes

## API design

R lists are represented by typed derived structures. Vector-valued prices are
returned in `option_prices`, optimizer information in `optimizer_result`, and
model calibrations in dedicated fit types. Optional R arguments become Fortran
optional arguments.

R's `optim` call is replaced by a self-contained Nelder-Mead implementation.
The original objective penalties remain in place. Requested Hessians are
computed by centered finite differences at the fitted parameters.

## Numerical support

The translation includes normal and lognormal functions, beta density and CDF,
Gaussian elimination, quadratic least squares, simple linear regression,
implied-volatility bisection, and finite-difference Hessians. No external
statistics or linear-algebra library is required.

## Deliberate corrections and safeguards

1. `compute.implied.volatility` was called by upstream `MOE` with a negative
   volatility lower bound. Volatility is now constrained to be positive, and
   the upper bracket is expanded when necessary. Failed brackets return NaN.
2. `mln.am.objective` allocated one high-strike put value with
   `numeric(length(len.ind.1.puts))`, which always has length one, and used the
   number of high strikes as an offset where the number of low strikes was
   required. The Fortran routine evaluates each actual strike directly and
   therefore avoids both indexing defects and the upstream sorted-strike
   assumption.
3. `price.shimko.option` numerically integrates a lognormal payoff while holding
   the strike-specific volatility fixed. The Fortran version evaluates the
   analytically equivalent BSM formula, avoiding unnecessary quadrature error.
4. Tiny tail probabilities in `price_am_option` are guarded to avoid division
   by zero in conditional moments.

## Presentation changes

The original `MOE` routine creates a multi-page PDF and two CSV files. The
Fortran `moe` routine performs the numerical fits and returns results without
I/O side effects. Density curves and predicted prices remain directly
available through the public density and pricing functions.

## Naming

Dots in R names become underscores. Fortran source uses lower case,
`implicit none`, `dp = kind(1.0d0)`, free-form source, and one statement per
line.
