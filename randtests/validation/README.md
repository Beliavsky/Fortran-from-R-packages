# Validation

From the package root:

```sh
Rscript validation/generate_reference.R
fpm build --profile debug
fpm test --profile debug
fpm run --example randomness_demo --profile debug
fpm build --profile release
fpm test --profile release
fpm run --example randomness_demo --profile release
```

R is used only for regenerating checked-in test fixtures. FPM tests and examples
run without R. `generate_reference.R` sources the unmodified pinned R snapshot;
it creates `r_reference.csv` and `test/r_reference.f90` reproducibly. Fixtures
contain 72 cases: eight method variants, three alternatives, and three sequences
(untied, tied, and decreasing). Each checks statistic, p-value, and sample size.
The p-value fixture explicitly caps upstream values at one; raw upstream values
remain visible in the CSV. Normal/Beta and unaffected exact runs computations
are checked at relative/absolute tolerance 3e-12.

Independent mathematical checks cover every binary sequence for all class
splits through total length 10. Those reference PMFs/CDFs validate runs mass,
CDF, upper and log tails, noninteger cutoffs, and quantile inequalities. Exact
Bartels histograms for n=1..7 are compared to a separate materialized permutation
algorithm with a squared-difference callback; exact n=9 verifies 9!=362880 total
permutations. Additional explicit regression cases cover the corrected n=1/n=2
Bartels masses, monotone n=3 exact test, inclusive right tail, capped Cox p-value,
and degenerate exact runs.

Validation also checks NaNs/infinities, insufficient/constant/tied samples,
invalid alternatives and methods, count overflow/resource bounds, invalid
quantile probabilities, empty/partial permutations, null callbacks, helper
statistics (including large-offset translation invariance), quantile CDF boundary
round trips, RNG support, and same-runtime fixed-seed reproducibility. No RNG
sequence equality with R or across different Fortran compilers is asserted.

The deterministic suite currently reports 20774 assertions. Runtime environment
used: R 4.5.0, gfortran 14.2.0, FPM 0.13.0, Linux. Windows is not available here;
no OS-specific source or linker dependency is required by this package.
