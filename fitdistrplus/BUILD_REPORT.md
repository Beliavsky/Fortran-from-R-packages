# Build report

- Translation version: 0.1.0
- Upstream version: 1.2-6
- Compiler: GNU Fortran (Debian 14.2.0-19) 14.2.0
- Fortran standard mode: Fortran 2018
- Checked suite: five of five test programs passed
- Optimized suite: five of five test programs passed
- Demonstration: passed; recovered Weibull shape 1.81290 and scale 2.49936 from seeded quantile data with true values 1.8 and 2.5
- Executable stack audit: GNU_STACK is RW, not executable
- FPM executable: unavailable
- FPM manifest: parsed successfully as TOML

The optimized build disables only GNU Fortran's optimizer-only equivalent of `-Wuninitialized` and `-Wmaybe-uninitialized` for allocatable
descriptors in local derived-type results. The checked build keeps all warnings
enabled and promoted to errors.
