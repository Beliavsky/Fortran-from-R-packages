# Porting notes

## Representation

The upstream package represents transforms as named R lists with attributes.
The Fortran port uses derived types such as `wavelet_transform`,
`packet_transform`, `wavelet_transform_2d`, and `wavelet_transform_3d`.
Coefficient arrays retain the same scale ordering as the R routines.

## Numerical implementation

- Native C DWT, inverse DWT, 2-D/3-D transform, and Hosking recursions were
  translated to free-form Fortran.
- The FFT is a self-contained radix-2 complex transform.
- DPSS tapers use a self-contained symmetric eigensolver.
- Long-memory likelihood routines use deterministic numerical searches rather
  than R optimizers.
- Random simulations can be seeded through `seed_rng` or routine-specific seed
  arguments. They do not reproduce R's Mersenne-Twister stream bit for bit.

## Boundary rules

Periodic and reflection boundaries are supported for the principal 1-D, 2-D,
and 3-D transforms. Two-dimensional wavelet-packet transforms currently use
periodic boundaries, matching the native packet decomposition implemented in
this port.

## Source parity and corrections

- `wavelet_filter` applies cascade letters from right to left, matching the
  upstream R implementation.
- `squared_gain` evaluates the complete cascaded filter rather than treating
  each filter letter independently.
- `cplxdual2D` and `icplxdual2D` in the supplied R source call `pm` with one
  argument although `pm` requires two. The port provides the tested
  `dualtree_2d`/`idualtree_2d` implementation instead of reproducing that
  broken call path.
- R plotting and list/name manipulation are intentionally not translated.

## Scope limits

The separate Hilbert wavelet-packet routine and its packet phase-shift helper
are not exposed as independent Fortran APIs. Equivalent ordinary packet,
Hilbert DWT/MODWT, and dual-tree analyses are available.
