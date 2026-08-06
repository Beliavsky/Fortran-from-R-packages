# Coverage

## Included

All computational routines exported by the package are represented:

- record counting;
- permutation averaging and confidence bounds;
- calibrated record-to-SNR splines;
- sample-size and tail corrections;
- normality testing;
- automatic tail-exponent fitting;
- final SNR and confidence estimates.

The three serialized R spline closures were decoded into their knot, value,
and cubic-coefficient arrays. Their IEEE binary64 values are embedded exactly.

## Replaced by typed Fortran structures

- R lists are represented by `r0_result` and `snr_result`;
- R missing-value filtering is implemented with IEEE finite checks;
- R/C++ random shuffling is represented by a local reproducible generator and
  Fisher-Yates shuffle.

## Omitted

- Rcpp registration and `.Call` wrappers;
- R package loading and unloading hooks;
- R help/database machinery;
- R serialized closure execution at runtime;
- bundled example datasets as compiled library inputs.

There is no plotting code in the upstream package.
