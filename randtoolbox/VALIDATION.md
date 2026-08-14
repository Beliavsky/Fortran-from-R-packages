# Validation

Validation was performed with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O2
```

Results:

- SFMT: all ten supported Mersenne exponents matched the independently compiled
  upstream C implementation for reference raw 32-bit output vectors.
- MT19937: 10,000 raw 32-bit outputs matched upstream C exactly: four scalar
  seeds (0, 1, 123456789, 0xffffffff), 2,000 outputs each, plus 2,000 outputs
  from the standard four-word array seed.
- Knuth TAOCP: 8,004 generated doubles matched upstream C exactly for seeds
  0, 1, 1302, and 123456789.
- Sobol: 23,120 32-bit direction-number words matched the bundled C table
  exactly for dimensions 1, 2, 3, 7, 32, and 1111 using 20 direction bits.
- General congruential RNG: 4,000 transitions matched an independent Python
  arbitrary-precision implementation, including non-power-of-two, power-of-two,
  and modulo-2^64 cases.
- WELL: the vendored implementation was previously validated bit-for-bit over
  272,000 raw outputs covering all 17 variants, four seeds, and 4,000 draws.
- Permanent FPM tests cover fixed SFMT/MT/WELL vectors, Knuth reproducibility,
  known Halton/Torus/Sobol points, bit conversions, Stirling numbers,
  permutation shape, collision counting, and chi-square frequency counts.

`fpm` itself was not installed in the validation environment.  The complete FPM
layout and `fpm.toml` were compiled directly with gfortran using the flags above.
