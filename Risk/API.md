# API reference

All public functionality is available through:

```fortran
use risk
```

All real arguments and results use `real(dp)`, where:

```fortran
integer, parameter :: dp = kind(1.0d0)
```

For unbounded support, pass `-huge(1.0_dp)` or `huge(1.0_dp)`.

## Distribution types

```fortran
type(normal_distribution)      :: d1  ! mu, sigma
type(lognormal_distribution)   :: d2  ! meanlog, sdlog
type(uniform_distribution)     :: d3  ! a, b
type(exponential_distribution) :: d4  ! rate
type(logistic_distribution)    :: d5  ! mu, scale
type(student_t_distribution)   :: d6  ! nu, mu, scale
```

Each extends `continuous_distribution` and implements:

```fortran
value = dist%pdf(x)
value = dist%cdf(x)
value = dist%quantile(p)
value = dist%lower_bound()
value = dist%upper_bound()
```

### Callback distribution

```fortran
type(callback_distribution) :: dist
call dist%initialize(pdf_proc, cdf_proc, quantile_proc, lower, upper)
```

Each callback has the interface:

```fortran
function callback(x) result(y)
   use risk_kinds, only : dp
   real(dp), intent(in) :: x
   real(dp) :: y
end function callback
```

## Risk measures

The argument `dist` below is `class(continuous_distribution)`.

```fortran
x = varg(dist, alpha)
x = esg(dist, alpha)
x = tcm(dist, alpha)
x = expp(dist, alpha, a, b)
x = bvar(dist, alpha, a)
x = epsg(dist, alpha)
x = expect(dist, a, b)
x = expvar(dist, alpha, a, b)
x = omegag(dist, alpha, a, b)
x = sortinog(dist, alpha, a, b)
x = kappag(dist, alpha, n, a, b)
x = wangg1(dist, alpha, a, b)
x = wangg2(dist, alpha, a, b)
x = stoneg1(dist, x0, k, a, b)
x = stoneg2(dist, x0, k, a, b)
x = luceg1(dist, a, b, aa, bb)
x = luceg2(dist, a, b, aa, bb)
x = luceg3(dist, a, b, aa, bb)
x = luceg4(dist, a, b, aa, bb)
x = saring1(dist, a, b, k, c)
x = saring2(dist, a, b, aa, bb1, bb2)
x = saring3(dist, a, b, aa, bb1, bb2)
x = bkg1(dist, alpha, a, b)
x = bkg2(dist, alpha, a, b)
x = bkg3(dist, alpha, a, b, beta)
x = bkg4(dist, alpha, a, b, beta)
```

For the following generic names, `alpha` may be either a scalar or a rank-one
array, and the result has the corresponding rank:

```text
varg esg tcm expp bvar epsg expvar omegag sortinog kappag
wangg1 wangg2 bkg1 bkg2 bkg3 bkg4
```

## Error results

Invalid domains return an IEEE quiet NaN. Examples include:

- probabilities outside the permitted range
- nonpositive scale parameters
- nonpositive Kappa order
- `saring3` with `aa = 1`, for which the source formula is singular
- callback distributions with missing procedure pointers
