# Upstream provenance

Translated from the supplied source archive of R package **GPareto 1.1.9**.
The upstream package is Copyright its listed authors/contributors and licensed
under GPL-3.

Important upstream native references retained under `upstream/native/`:

- `domination.cpp`: Kung-style non-domination and matrix distance kernels.
- `EHI_2d_rcpp.cpp`: analytical two-objective expected hypervolume improvement.

GPareto depends on DiceKriging. The required DiceKriging computational modules
are vendored from `DiceKriging-fortran-v0.1.0`, electing DiceKriging's GPL-3
license option so the combined GPareto distribution is GPL-3. See
`licenses/DiceKriging-GPL-3.0.txt`.

The original GPareto `DESCRIPTION` and `NAMESPACE` are retained in `upstream/`.
