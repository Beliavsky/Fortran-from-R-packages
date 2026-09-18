# kza 4.2.0

This release fixes several long-standing bugs found in a source audit of
the 4.1.x code, adds a regression test suite, and modernizes the package
metadata. **Several fixes change numerical results**, most significantly
for two-dimensional (matrix) input. If you have published or archived
analyses that used `kza()` on matrices or 3-D arrays, or `kzsv()`, we
recommend re-running them under this version; 1-D (vector and time
series) results are unchanged except where noted.

## Results-changing bug fixes

* **`kza()`'s `m` now means the full window width, matching `kz()`.**
  Previously `kza()` used `m` itself as the per-side radius while `kz()`
  used `floor(m/2)`, so the adaptive filter's window was roughly twice
  its own baseline's for the same `m` — inconsistent with the published
  algorithm, which drives the baseline and the adaptive offsets with one
  half-width. `kza()` and `kzsv()` now use `floor(m/2)` as the radius.
  Smoothing for a given `m` is therefore about half as wide as before;
  double `m` to reproduce old results.

* **2-D `kza()` now actually iterates.** The 2-D code path never carried
  one iteration's result into the next, so matrix input silently
  returned the `k = 1` result regardless of `k` (the 1-D and 3-D paths
  iterated correctly). Matrix results with `k > 1` will differ — they
  are now smoother in homogeneous regions, as documented.

* **2-D `kza()` addressed the matrix transposed inside its box
  average.** On square matrices this self-consistently averaged the
  transposed pixel's neighborhood, shifting results modestly; on
  non-square matrices it read far outside the matrix, so output could
  contain arbitrary values that varied from run to run. Non-square 2-D
  results from 4.1.x and earlier should be considered undefined; square
  2-D results will shift slightly under this version.

* **2-D `kza()` no longer shrinks the vertical window on tall, narrow
  matrices.** A leftover bounds check clamped the row-direction window
  against the number of *columns*, wrongly limiting smoothing whenever
  the window exceeded `ncol(x)`. The filter now treats rows and columns
  symmetrically (for a scalar window, filtering commutes with matrix
  transposition).

* **`kza()` on matrices and arrays no longer overwrites cells with
  `NA`.** The `impute_tails = FALSE` default marked the first and last
  `m` elements `NA` by linear index regardless of dimensionality, which
  blanked meaningless column-major positions in 2-D/3-D output. Tail
  marking now applies to 1-D input only; 1-D behavior is unchanged.

* **`kzsv()` now honors `min_size` and ignores `k`.** The arguments
  passed to the compiled code were shifted by one position: the
  iteration count was consumed as the minimum window length and the
  fitted `min_size` was never passed at all. Sample-variance results
  will differ accordingly.

* **3-D `kza()` no longer skips averaging when a window collapses to
  width 1 in some dimension.** Affected voxels previously received the
  raw center value instead of the average over the rest of the box; in
  the extreme case, a single-slice `(n, m, 1)` array was returned
  completely unfiltered. Results change at faces and near strong
  gradients, and single-slice arrays are now genuinely filtered.

## New capabilities (all opt-in unless noted)

* `kza(..., symmetrize = TRUE)` (matrix input): averages the filter over
  the four 90-degree rotations of the input. The adaptive head/tail rule
  can misallocate the window immediately beside a sharp, high-contrast,
  axis-aligned edge, so the plain filter does not commute with 90-degree
  rotation; rotation averaging cancels that row/column bias -- exactly,
  for square input -- at four times the compute.

* `kza(..., normalize = "quantile")`: scales the adaptive shrink factor
  against the 99th percentile of the difference metric instead of its
  literal maximum, so a handful of extreme pixels do not single-handedly
  set the smoothing scale for the whole image. Falls back to the maximum
  automatically when the quantile is zero (sparse structure on an
  otherwise flat field would otherwise be smoothed away). The default
  (`"max"`) is byte-identical to previous behavior.

* `rlv()` reimplemented on cumulative-sum tables: ~35x faster on a
  256x256 image at `krnl = 9`. **The default boundary policy changes**:
  `pad = "clamp"` computes edge windows from the cells actually present,
  where the old code conceptually zero-padded the data and mixed zeros
  into border variances; interior cells are unchanged, and
  `pad = "zero"` reproduces the old behavior exactly. An even `krnl`
  greater than 1 -- which the old code silently mis-indexed -- is now an
  error; `krnl = 1` keeps the one-sided max-of-corners mode.

## Other fixes

* The adaptive shrink factor is guarded against a zero difference-metric
  maximum: constant input previously computed `0/0` and relied on
  undefined C behavior downstream; it now passes through exactly.

* `kzsv()` accepts objects that subclass `"kza"` (`inherits()` instead
  of `class() ==`), and the same fix applies to the `kzp` check in the
  periodogram code.

* The native-routine registration for `R_kzsv` pointed at the wrong
  function with the wrong argument count (latent — no R code called that
  entry point).

* The 1996 reference in `?kza` links the paper's DOI; the old AMS search
  URL had gone dead.

* Three-dimensional `kza()` no longer reads outside its internal
  working arrays when computing the adaptive window near the far y-edge
  of any array with more than one row, or anywhere in a single-slice,
  single-row, or single-column 3-D array. Found via AddressSanitizer;
  like the adaptive-shrink guard above, this was undefined behavior that
  could bias the adaptive window's orientation at those positions rather
  than corrupt results outright, so most affected results are unchanged
  or shift only slightly at those edges.

## Infrastructure

* New regression test suite (testthat, 35 tests). Every bug fix above
  was first demonstrated to fail against an unmodified 4.1.0.1 build,
  except the two most recent undefined-behavior fixes (the adaptive
  shrink guard and the kza3d edge-array reads), which have no
  deterministic red/green diff on a plain build; those were confirmed
  against the baseline with AddressSanitizer instead.

* `Authors@R` metadata and `.Rbuildignore` added.
