# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of zigg 0.0.2.
- Ported normal, exponential, and uniform Ziggurat generators.
- Ported 32-bit KISS/SHR3 state engine with explicit modulo-2^32 arithmetic.
- Added module-level API compatible with the four R exports.
- Added reusable `ziggurat_rng` derived type and state save/restore.
- Added exact C++ stream regression tests and statistical smoke tests.
