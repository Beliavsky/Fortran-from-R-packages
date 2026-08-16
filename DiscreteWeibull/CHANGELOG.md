# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of `DiscreteWeibull` 1.1.
- Ported complete type-I and type-III distribution APIs.
- Ported moment, likelihood, loss, and all three estimation methods.
- Integrated supplied `Rsolnp-fortran` for type-I moment fitting.
- Added observed-information covariance.
- Corrected the documented type-III hazard formula's missing multiplication by `c`.
