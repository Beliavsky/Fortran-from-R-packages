# Licensing and attribution

This project is a translation of computational code from R package `nnet`
version 7.3-21.

Upstream authors/copyright holders:

- W. N. Venables
- B. D. Ripley
- portions distributed with the first edition were also copyright Springer-
  Verlag New York Inc., with rights assigned as described in upstream
  `LICENCE.note`.

The upstream package is licensed, at the recipient's option, under GNU GPL
version 2 or GNU GPL version 3.  The translated `nnet`-derived Fortran files
carry `GPL-2.0-only OR GPL-3.0-only` identifiers.  Copies of GPL-2 and GPL-3 are
included under `licenses/`.

The supplied `r_mod.f90` remains under the MIT license.  Its original is
retained under `upstream/` and its MIT notice is included separately.

The multinomial Fisher-information R code notes that the analytic information
calculation was contributed by David Firth and cites the formulation in
T. Amemiya, *Advanced Econometrics* (1985), pp. 295-296.

The neural-network Hessian code cites B. D. Ripley (1996), *Pattern Recognition
and Neural Networks*, p. 152.  The package citation in `upstream/inst/CITATION`
and original source headers are retained for complete provenance.
