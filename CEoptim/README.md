# CEoptim-fortran

Modern Fortran/FPM translation of the computational core of the R package
**CEoptim 1.3**, a Cross-Entropy optimization package.

## Implemented functionality

The library translates the numerical code in `R/CEoptim.R`, `R/rtmvnorm.R`,
and `R/dirichletrnd.R`:

- continuous Cross-Entropy optimization with independent Gaussian sampling;
- categorical/discrete Cross-Entropy optimization;
- mixed continuous/categorical optimization;
- minimization and maximization;
- elite selection using `rho`;
- separate smoothing of means, standard deviations, and categorical
  probabilities;
- convergence thresholds for standard deviations and categorical
  probabilities;
- the package's no-improvement termination rule;
- linear continuous constraints `A x <= b`;
- truncated multivariate-normal sampling by accept/reject with adaptive batch
  sizing and Gibbs fallback;
- multivariate-normal sampling from positive-semidefinite covariance matrices;
- univariate truncated-normal sampling;
- Dirichlet random variates using gamma normalization;
- CE iteration history and categorical probability history;
- package-compatible `nfe = niter*N`, plus `actual_nfe`, which also counts the
  initial population evaluation.

R S3 printing, R list/formula handling, imported graph plotting, R's global RNG,
and package data objects are not computational algorithms and are not ported.
The user's objective callback can capture fixed parameters directly, which is
the Fortran equivalent of R's `f.arg` mechanism.

## Public API

The umbrella module is `ceoptim`.

Main types:

```fortran
use ceoptim

type(ce_control)             :: control
type(ce_continuous_control)  :: continuous
type(ce_discrete_control)    :: discrete
type(ce_result)              :: result
```

The optimizer is called as

```fortran
call ce_optimize(objective, result, control, continuous, discrete)
```

Either `continuous` or `discrete` may be omitted.

The objective has the interface

```fortran
function objective(xc, xd) result(value)
   use ceoptim, only : dp
   real(dp), intent(in) :: xc(:)
   integer, intent(in)  :: xd(:)
   real(dp) :: value
end function objective
```

Categorical values follow CEoptim exactly: a variable with `k` categories takes
integer values `0,1,...,k-1`.

For nonuniform initial categorical probabilities, allocate
`discrete%probs(maxval(categories), q)` and put each variable's probabilities
in the first `categories(j)` rows of column `j`.

## Linear constraints

For continuous variables,

```fortran
continuous%con_mat = A
continuous%con_vec = b
```

represents

```text
A x <= b
```

The translated `rtmvnorm` first tries untruncated multivariate-normal proposals,
as CEoptim does. If the empirical acceptance rate becomes too small, it switches
to component-wise Gibbs sampling. The Gibbs implementation uses the equivalent
precision-matrix conditional-normal formulas rather than reconstructing every
leave-one-coordinate covariance inverse separately.

## Random-number generation

The original R package uses R's RNG, `MASS::mvrnorm`, `msm::rtnorm`, and
`stats::rgamma`. This port uses a self-contained Park-Miller uniform generator,
Box-Muller normals, Marsaglia-Tsang gamma sampling, and a robust truncated-normal
sampler. Set `control%seed` or call `rng_seed` for reproducibility.

A numerical seed therefore does **not** reproduce R's exact random stream, but
it does reproduce a Fortran run of this library.

## Build with FPM

```text
fpm build
fpm test
fpm run --example peaks
fpm run --example mixed
```

BLAS and LAPACK are linked through `fpm.toml`.

FPM was not installed in the translation environment, so the FPM source layout
was validated with equivalent direct GNU Fortran builds against system
BLAS/LAPACK.

## Validation

The tests include:

1. the package documentation's two-dimensional Peaks maximization problem;
2. exact categorical optimization;
3. mixed continuous/discrete optimization under `A x <= b` constraints;
4. Dirichlet sampling;
5. both accept/reject and forced-Gibbs truncated multivariate-normal paths.

With GNU Fortran 14.2.0, the Peaks regression converges to approximately

```text
maximum = 8.10621359
x       = (-0.00932, 1.58138)
```

which is the known global peak of the standard Peaks test function to the
accuracy expected from stochastic CE optimization.

## Source layout

```text
src/        Fortran library
test/       regression tests
example/    FPM examples
original/   retained CEoptim package metadata and R sources
COPYING     GPL v2 text
```

See `TRANSLATION_NOTES.md` for algorithm-by-algorithm correspondence and known
behavioral differences.
