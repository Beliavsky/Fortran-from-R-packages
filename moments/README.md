# moments-fortran

Modern Fortran 2018 translation of the computational code in the R package
`moments` 0.14.1.

The library computes raw, central, and absolute moments; converts complete
moment sequences; evaluates cumulants; calculates skewness and Pearson/Geary
kurtosis; and implements the D'Agostino, Anscombe-Glynn, Bonett-Seier, and
Jarque-Bera tests.

## Build with FPM

```text
fpm build
fpm test
fpm run
```

The public module is:

```fortran
use moments
```

Moment sequences use normal Fortran indexing: array element `k + 1` contains
the order-`k` moment or cumulant. Vector and matrix overloads are available;
matrix calculations operate column by column.

## Example

```fortran
use moments, only : dp, all_moments, skewness, jarque_test, moments_test_result
real(dp) :: x(5)
real(dp), allocatable :: mu(:)
type(moments_test_result) :: jb

x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 8.0_dp]
mu = all_moments(x, order_max=4)
jb = jarque_test(x)
print *, mu(2), skewness(x), jb%p_value
```

## Scope

All 12 namespace exports are represented. R vectors, matrices, and data frames
become explicit real arrays. R `htest` objects become `moments_test_result`.
There is no plotting or external numerical dependency.

The upstream package is licensed under GPL version 2 or later. This translation
is distributed under `GPL-2.0-or-later`.
