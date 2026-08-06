# Changes in Version 1.4-0 (DA)
- span_as(): new control parameters `B` and `seed`, de-randomising the L > 0
  test. For L > 0 the score is multiplied by a random weight vector that is
  SHARED by every test asset, so its effect does not average out across the
  cross-section: a single draw is a single realisation of the test. The test is
  valid for any fixed draw (its size is correct), but the p-value returned on one
  data set can move by orders of magnitude from draw to draw -- most visibly when
  T is short, where floor(T^(1/3)) leaves few subseries. `B > 1` draws B
  independent weight vectors and merges the per-asset p-values with the Cauchy
  rule, which is valid under arbitrary dependence and is not carried by a single
  extreme draw. `B = 1` (the default) reproduces the previous behaviour exactly;
  `B = 100` or more is recommended for any reported application. The seed of the
  first draw is now user-visible (`seed`, default 123, draw b uses seed + b - 1)
  instead of hard-coded. Consecutive seeds were kept, rather than an explicit
  substream generator, so that `B = 1` reproduces the historical result exactly;
  the draws they produce are empirically indistinguishable from independent
  (mean absolute pairwise correlation 0.051 over 100 seeds at T = 250, against
  0.050 expected), and this is documented under "Choosing B".
- span_as(): `B` and `seed` must now be whole numbers. Previously a fractional
  `B = 1.9` passed validation and was silently truncated to a single draw.
- span_gl_a() / span_gl_ad(): a singular benchmark design (e.g. collinear
  benchmarks) now returns NA p-values instead of raising, matching the
  convention already used by span_grs(), span_f2() and span_py().
- The Cauchy combination now returns NA, not a silent NaN, when it is handed no
  usable p-values (every test asset missing).
- The subseries t-test now fails with an explicit message when the sample is too
  short to form two subseries, instead of dividing by zero degrees of freedom.

# Changes in Version 1.3-0 (DA)
- span_gl_a() / span_gl_ad(): the C++ sign-flip kernel now forms the restricted
  residuals via Ehat0s = Ehat1 + (XX %*% premult) %*% (H %*% Bhat1 - C), avoiding
  the per-simulation T x (K+1) x N restricted-least-squares matmul. This removes
  one of the three large per-simulation matrix products, cutting the GL kernel
  time by roughly 20% (about 1.2x at K = 100, N = 200, T = 250), with
  bit-identical p-values.
- span_gl_a() and span_gl_ad() now default to do_trace = FALSE (no console
  output unless requested).
- span_grs(), span_f2() and span_py() now return NA (like the other tests)
  instead of raising an error when a covariance/design matrix is singular
  (e.g. collinear benchmarks).
- span_simulate() validates its arguments (dimensions, correlation ranges,
  degrees of freedom, GARCH stationarity).
- span_as() is faster: the internal subseries t-test (f_ttest) no longer calls
  stats::t.test() once per column (whose input-checking dominated the runtime).
  The one-sample two-sided p-values are now computed directly and vectorised
  across columns, giving identical results (to floating-point) at about 2.5x
  lower cost for span_as().
- NEW: span_simulate() generates benchmark and test-asset returns for spanning
  size/power studies from a factor model with a controllable spanning violation
  (ncp). Innovations may be normal, Student-t, or skew-t, with iid, AR(1),
  GARCH(1,1), or AR-GARCH dynamics; the `dgp` argument selects the twelve
  processes of Ardia and Sessinou (2025). It is a fast, validated drop-in for
  ad-hoc fGarch::garchSim() loops (roughly 30-45x faster for the GARCH DGPs) and
  matches the reference processes in distribution, dynamics, and test size/power.
  The skew-t innovations use a base-R re-implementation of the Fernandez-Steel
  standardised skew-t (identical draws to fGarch::rsstd given the same RNG
  state), so the package has no external simulation dependency. Its GARCH(1,1) recursion runs
  in C++ (bit-for-bit identical to the R loop; floating-point contraction is
  disabled so the fused multiply-add does not alter the roundings), which about
  halves the data-generation time for the GARCH DGPs.
- span_gl_a() and span_gl_ad() now run their sign-flip Monte Carlo simulations
  in C++ (via Rcpp / RcppArmadillo). The random signs and the tie-breaking
  uniforms are still drawn in R, so the RNG stream -- and hence every p-value and
  decision -- is unchanged (verified identical to 1.2-1 across dimensions,
  totsim, and seeds); only the per-simulation SSR / F_max computation moved to a
  streaming C++ kernel that never forms the large T x (N*totsim) intermediates.
  About 3x faster than the 1.2-1 R implementation at large K and N, with much
  lower memory use, which also improves the parallel scaling of simulation
  studies. Installing from source now requires a C++ compiler (Rcpp and
  RcppArmadillo were added as build dependencies).

# Changes in Version 1.2-1 (DA)
- span_as() is now vectorized across the test cross-section: benchmark-only QR
  decompositions are formed once and the swap regressions are obtained by
  Frisch-Waugh partialling, replacing the per-asset loop of full QR
  factorizations. Output is numerically identical to 1.2-0 (verified); span_as
  is roughly 12-15x faster for moderate-to-large test sets.
- span_gl_a() and span_gl_ad() streamlined. The sign-flip simulations are now
  assembled directly as T x (N*totsim) matrices (recycling the residuals and
  fitted values across simulations, gathering the sign columns) instead of
  building a T x N x totsim array and transposing it with aperm(), which removed
  the dominant cost at large N. In addition the balanced-MC restricted SSR is
  constant across sign-flips (it equals the raw restricted SSR, since squaring
  removes the flipped sign) so it is computed once; redundant array/matrix
  reshapes in the constrained-estimate step were collapsed to matrix products;
  and unused Sigma allocations were removed. Output is bit-for-bit identical to
  1.2-0 (verified via identical()); GL is roughly 1.5-1.7x faster, most visibly
  for large test cross-sections (e.g. K = 100, N = 200).

# Changes in Version 1.2-0 (DA)
- NEW: span_as(), the Ardia-Sessinou subseries-based Cauchy Combination Test (CCT)
  for high-dimensional spanning, is now available and exported. It returns combined
  p-values for the alpha (CCTa), variance/slope (CCTd), and joint (CCTad) nulls,
  robust to serial/cross-sectional dependence and conditional heteroskedasticity.
  Internal helper f_getpv() implements the per-asset computation.
- span_bj now implements the Britten-Jones (1999) tangency-spanning test correctly
  (regression of ones on raw returns); the previous version had no power against
  alpha alternatives
- span_py now returns NA instead of erroring when there is a single test asset (N = 1)
- f_prods no longer overwrites the caller's global RNG state (.Random.seed is
  saved and restored)
- span_gl_ad trace now prints the decision label instead of the numeric code
- tail p-values in span_grs and span_py use lower.tail = FALSE for better precision
- added @importFrom Rdpack reprompt to silence the R CMD check NOTE on unused Imports
- documentation: benchmark/test dimensions in span_km, span_gl_a, and span_gl_ad
  aligned with the rest of the package (R1 is T x K, R2 is T x N)
- added correctness (size/power) tests for span_bj and an N = 1 test for span_py
- added tests for span_as and f_getpv

# Changes in Version 1.1-3 (DA)
- seed is finally removed from GL functions

# Changes in Version 1.1-2 (BS)
- DESCRIPTION modified according to CRAN guidelines
- seed is now an input for GL functions
- error fixed in the test file of span_bj

# Changes in Version 1.1-1 (BS)
- Function span_km and span_bj were testing the same

# Changes in Version 1.1-0 (DA)
- First release public version
