# Translation coverage

The upstream NAMESPACE contains 272 exports.

Covered computational areas:

- 239 compiled special-function/vector entry points from the upstream C layer
- Airy functions and zeros
- cylindrical, modified and spherical Bessel families and zeros/sequences
- Clausen, Dawson, Debye and dilogarithms
- Coulomb wave functions and coupling coefficients
- complete/incomplete Carlson and Legendre elliptic integrals and Jacobi functions
- error functions and hazard
- exponential integrals, sine/cosine/hyperbolic integrals and atan integral
- Fermi-Dirac integrals
- gamma/beta/incomplete gamma, factorial/binomial and Pochhammer functions
- Gegenbauer and Laguerre polynomials
- hypergeometric functions
- Lambert W
- Legendre and associated/spherical/conical/H3d functions
- logarithmic and complex logarithmic functions
- integer powers
- psi/polygamma functions
- synchrotron and transport functions
- trigonometric and complex trigonometric functions
- zeta/eta/Hurwitz zeta families
- polynomial evaluation
- all 14 upstream RNG types
- Sobol and Niederreiter-2 QRNGs
- complex Jacobi aliases built in R

Not reproduced as active computations:

- R attribute preservation and automatic argument recycling (`process.args`)
- R package load hooks
- R external-pointer classes
- `multimin*`, because the upstream 2.1-9 R functions themselves are disabled and immediately stop

The project intentionally links to GNU GSL rather than copying or independently rewriting the GSL numerical algorithms.
