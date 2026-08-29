# Validation

The port is validated by five retained test programs:

1. `test_core` -- independent density/CDF/mean/variance values, untruncated
   identities, support behavior and the package's extreme-tail fallback.
2. `test_quantile` -- `p(q(p))` inversion over finite, left-truncated,
   right-truncated and shifted/scaled cases.
3. `test_random` -- Monte Carlo mean/variance and support for a nontrivial
   two-sided distribution.
4. `test_recycle` -- R-style recycling and vector generic interfaces.
5. `test_sampler_branches` -- Monte Carlo checks exercising ordinary Normal,
   central two-sided, half-normal, uniform/exponential and symmetry branches.

Independent reference values in `test_core` were generated from the standard
truncated-normal formulas using SciPy's Normal CDF/PDF, not from the translated
Fortran implementation.

The strict validation build uses GNU Fortran with Fortran 2018,
`-fcheck=all`, and `-Werror=implicit-interface` and links BLAS/LAPACK for the
supplied `r_mod` object.
