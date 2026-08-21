# BGFD-fortran

Modern Fortran translation of the computational code in the R package **BGFD
0.1 (Bell-G and Complementary Bell-G Family of Distributions)**.

The project is an FPM library, uses `dp = kind(1.0d0)`, requires no R runtime,
and preserves the upstream GPL (>= 2) licensing as GPL-2.0-or-later.

## Implemented families

Eight baseline distributions are available under both Bell-G and complementary
Bell-G transforms:

- exponential (`bell_e`, `cbell_e`)
- exponentiated exponential (`bell_ee`, `cbell_ee`)
- Weibull (`bell_w`, `cbell_w`)
- exponentiated Weibull (`bell_ew`, `cbell_ew`)
- Fisk/log-logistic (`bell_f`, `cbell_f`)
- Lomax (`bell_l`, `cbell_l`)
- Burr XII (`bell_b`, `cbell_b`)
- Burr X (`bell_bx`, `cbell_bx`)

Each family has `d_*`, `p_*`, `q_*`, `r_*`, `s_*`, `h_*`, and `m_*` routines.
The generic `bgfd_pdf`, `bgfd_cdf`, `bgfd_quantile`, `bgfd_random`,
`bgfd_survival`, and `bgfd_hazard` routines accept run-time family IDs.

## Example

```fortran
program demo
    use bgfd
    implicit none
    real(dp) :: x

    x = 1.3_dp
    print *, d_bell_w(x, 0.9_dp, 1.6_dp, 0.7_dp)
    print *, p_bell_w(x, 0.9_dp, 1.6_dp, 0.7_dp)
    print *, q_bell_w(0.37_dp, 0.9_dp, 1.6_dp, 0.7_dp)
end program demo
```

For MLE fitting:

```fortran
use bgfd
real(dp) :: x(100)
type(bgfd_fit_result) :: fit

! fill x ...
call m_bell_e(x, 1.0_dp, 1.0_dp, fit)
print *, fit%params
print *, fit%se
print *, fit%aic, fit%bic, fit%ks_statistic, fit%ks_pvalue
```

`bgfd_fit_result` also contains AICc, HQIC, Cramer-von Mises,
Anderson-Darling, `-2 log L`, convergence status, and density/CDF endpoint
checks.

## Build

With FPM:

```text
fpm build
fpm test
fpm run --example example_bgfd
```

A strict GNU Fortran validation script is also included:

```text
./run_tests.sh
```

It uses:

```text
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all
```

## Numerical corrections

The R package's four Weibull-based quantile routines use `beta` where the
inverse CDF requires `1/beta`. This is corrected here. The R survival/hazard log
flags also combine log-probabilities with ordinary probabilities; the Fortran
API gives those flags conventional mathematical meanings. See
`PORTING_NOTES.md` for details.

## Attribution

The original BGFD authors are Michail Tsagris, Muhammad Imran, and M.H. Tahir.
The original supplied R package is retained under `upstream/`. The
GPL-compatible AdequacyModel Fortran translation used for MLE and goodness of
fit is vendored in the source tree with its notices under `third_party/`.
