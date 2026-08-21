# elliptic-fortran

Modern Fortran 2018 translation of the computational code in Robin K. S. Hankin's
R package **elliptic 1.5-1**, "Weierstrass and Jacobi Elliptic Functions".

The translation is self-contained and uses no external numerical library.  It is
organized as an FPM package and uses double-precision complex arithmetic throughout.

## Main coverage

- complete elliptic integral `k_complete()` by complex AGM, `nome()` and `nome_k()`;
- Jacobi theta functions theta1--theta4 and the first three derivatives of theta1;
- Neville theta functions and the twelve Jacobi elliptic quotients (`sn`, `cn`, `dn`,
  `ns`, `nc`, `nd`, `sc`, `sd`, `cs`, `cd`, `ds`, `dc`);
- Abramowitz-Stegun identities/evaluation formulas 16.28, 16.36, 16.37, and 16.38;
- Weierstrass invariants, roots, primitive half periods, period reduction, and parameter
  conversion;
- Weierstrass `wp`, `wp_prime`, `weierstrass_sigma`, and `weierstrass_zeta`;
- Coqueraux duplication and Laurent-series implementations for P, P', sigma, sigma',
  and zeta;
- Eisenstein/invariant calculations (`g2`, `g3`) using theta, lattice, divisor, fixed,
  and Lambert-series forms;
- Dedekind eta, modular J and lambda;
- Farey/unimodular, Mobius-transform, lattice, congruence, prime-factorization,
  divisor, totient, Mobius, and Liouville utilities;
- complex contour integration, polygonal contour integration, residues, and complex
  Newton-Raphson.

R graphics (`view`, `latplot`) are intentionally omitted.  `P.pari()` is also omitted:
it is a shell wrapper around the external `pari/gp` executable, not an elliptic
algorithm implemented by the package.  R replacement/class helpers (`Re<-`, `Im<-`,
class attributes) are represented by ordinary Fortran complex arithmetic and the
`elliptic_parameters` derived type.

## Build

With FPM:

```text
fpm build
fpm test
fpm run --example basic
```

Or directly with GNU Fortran, compile the modules in this order:

```text
elliptic_kinds.f90
elliptic_theta.f90
elliptic_arithmetic.f90
elliptic_numeric.f90
elliptic_weierstrass.f90
elliptic_modular.f90
elliptic.f90
```

The release was validated with GNU Fortran 14.2 using Fortran 2018, `-O2`,
`-Wall -Wextra -Wimplicit-interface -Werror`, and `-fcheck=all`.

## License and upstream attribution

The upstream package is GPL-2.  This translation is distributed under
GPL-2.0-only.  See `LICENSE`, `UPSTREAM.md`, and `upstream/`.

Please cite the original work:

Robin K. S. Hankin (2006), "Introducing elliptic, an R package for elliptic and
modular functions", *Journal of Statistical Software* 15(7).
