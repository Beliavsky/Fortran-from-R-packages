# Porting notes

The translation preserves numerical algorithms and the MIT license. R-specific
S3/vctrs classes, tibbles, plotting, documentation-site files, and data-frame
adapters are outside the Fortran scope. Arrays use the Fortran convention
`x(n_scenarios,n_variables)`.

Equality-constrained entropy pooling is solved in the dual with damped Newton
iterations. A separate penalty/mirror-descent routine supports mixed linear
inequality and equality constraints without external optimization libraries.
