# rmutil-fortran

Modern Fortran 2018 translation of the computational portions of the R package
`rmutil` 1.1-10 (Bruce Swihart and Jim Lindsey).

This project is intended to expose the numerical/statistical algorithms in a
small FPM library.  R's S3 data objects, formula evaluator, printing, file I/O,
and graphics are not reproduced.

## Implemented computational surface

### Probability distributions

The complete exported `d/p/q/r` distribution surface in `rmutil::dist.r` is
represented:

- inverse Gaussian
- Laplace
- Levy
- Pareto
- simplex
- two-sided power
- Box-Cox
- Burr
- generalized extreme value
- generalized gamma
- generalized inverse Gaussian
- generalized logistic
- generalized Weibull
- Hjorth
- power exponential
- skew Laplace
- beta-binomial
- double binomial
- multiplicative binomial
- double Poisson
- multiplicative Poisson
- PVF Poisson
- gamma count
- Consul

Fortran names follow the R names (for example `pinvgauss`, `dggamma`,
`qdoublepois`, and `rconsul`). Random generators return allocatable vectors.
Optional density argument `log_p=.true.` corresponds to R's `log=TRUE`.

### Numerical methods

- `integrate_romberg`: one-dimensional Romberg integration, including infinite
  limits.
- `toms614_integrate`: direct modern Fortran port of ACM TOMS Algorithm 614
  (INTHP) carried by upstream `rmutil`.
- `integrate_adaptive`: adaptive Simpson integration.
- `integrate_2d`: mapped Gauss-Legendre product quadrature for finite or
  infinite rectangular domains.
- `runge_kutta`: scalar fourth-order Runge-Kutta solver.
- `matrix_exp`: self-contained scaling-and-squaring matrix exponential.
- `lin_diff_eqn`: autonomous linear ODE solution using the matrix exponential.
- `gauss_hermite`: Gauss-Hermite nodes and normalized weights.

### Repeated-measures/data numerical helpers

- `gettvc`: last-observation-carried-forward alignment of irregular
  time-varying covariates, including the upstream tie convention.
- `contrast_mean`: numerical equivalent of `contr.mean`.
- `group_sum` and `group_mean`: common numerical cases of `capply`.

### Pharmacokinetic/pharmacodynamic mean functions

All functions from upstream `R/pkpd.r` are included, with dots mapped to
underscores:

- `mu1_0o1c`, `mu1_1o1c`, `mu1_1o2c`, `mu1_1o2cl`, `mu1_1o2cc`
- `mu2_0o1c`, `mu2_0o2c1`, `mu2_0o2c2`, `mu2_1o1c`
- `mu2_0o1cfp`, `mu2_0o2c1fp`, `mu2_0o2c2fp`, `mu2_1o1cfp`

## Deliberately omitted R infrastructure

The following are principally R object-language/presentation facilities rather
than portable numerical algorithms and are not cloned:

- S3 `response`, `repeated`, `tccov`, `tvcov`, `formulafn`, profile, and `gnlm`
  method infrastructure
- `finterp`, `fnenvir`, `wr`, and R formula/model-frame evaluation
- `restovec`, `tcctomat`, `tvctomat`, `dftorep`, `rmna`, and `lvna` object/data
  constructors and class-preserving NA manipulation
- `read.*` functions
- all print and plot methods

The general real-valued matrix power operator `%^%` is also not reproduced;
`matrix_power_integer` is provided for integer powers. `matrix_exp` covers the
matrix exponential used by `lin.diff.eqn` without requiring LAPACK.

## Known upstream compatibility choices

Two upstream implementation details deserve explicit mention:

1. Upstream's multiplicative-binomial C implementation uses different
   interaction expressions in its density and cumulative paths. This port
   preserves those two paths rather than silently redefining the distribution.
2. The upstream PVF-Poisson C cumulative loop stops at `j < q`, so `P(Y<=0)`
   evaluates as zero although the documented routine is a CDF. This port uses
   the coherent inclusive cumulative definition `sum_{j=0}^q P(Y=j)` and
   defines the quantile from that corrected CDF.

The upstream skew-Laplace quantile uses `0.5` as its branch point even when the
asymmetry parameter makes `F(m) != 0.5`; the Fortran port preserves that
published implementation.

## Build

With FPM:

```text
fpm build
fpm test
fpm run --example
```

The sources have also been validated directly with GNU Fortran using Fortran
2018, explicit-interface errors, and runtime checking. See
`docs/VALIDATION.md`.

## Minimal use

```fortran
program demo
   use rmutil
   implicit none
   real(dp) :: p, ans

   p = pinvgauss(1.4_dp, 2.0_dp, 0.3_dp)
   ans = integrate_romberg(f, 0.0_dp, pi)
   print *, p, ans
contains
   function f(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = sin(x)
   end function f
end program demo
```

## License and provenance

The upstream package is GPL-2-or-later. This translation is distributed under
the same terms. `COPYING` contains GPLv2, and `upstream/` retains the original
package metadata plus the upstream TOMS614 and `gettvc` source files used for
provenance. TOMS Algorithm 614 is credited upstream to K. Sikorski, F. Stenger,
and J. Schwing.
