# Validation

The release is compiled with GNU Fortran using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

Regression programs cover:

1. linear total-least-squares solution against the analytic TLS line;
2. exact-line recovery and fixed-parameter behavior;
3. bounds, weights, and log-likelihood helpers;
4. the `window` path used when the number of observations exceeds 25;
5. the upstream Chwirut2 NIST data set.

The examples cover a noisy line and Chwirut2.

The translated library invokes user model callbacks only through explicit typed
module-level interfaces.  It does not rely on an internal procedure directly
calling a host-associated user procedure dummy, avoiding the portability issue
seen with some older Windows gfortran builds.
