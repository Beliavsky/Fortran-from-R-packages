# Porting notes

## Native C kernels

The original C formulas for circular cross-correlations, rank statistics,
bias corrections, Fisher combinations, cumulants, and Edgeworth
approximations are translated directly. The 1,001-point empirical CDF
tables from `cvmF50.c` and `cvmF100.c` are preserved as Fortran parameter
arrays. As upstream, samples below 61 observations use the 50-observation
table and larger samples use the 100-observation table.

## MixedIndTests dependency

The R package calls `MixedIndTests::EstDepMoebius` for Spearman, van der
Waerden, and Savage cross-dependence. The necessary score construction and
normal quantile code are adapted into this project, making the FPM library
self-contained.

## Lag representation

The C routines emit statistics in an internal traversal order and the R
wrappers rearrange them. The Fortran code computes directly in the final
natural order. Circular shifts are mathematically equivalent to the
original permutations.

## Corrected wrapper defect

In the negative-lag branch of upstream `crossdep_3series`, the `y`/`z`
pair uses a stale `y1` value from an earlier loop. The Fortran routine uses
the intended lag-specific circular shift for every pair. This affects only
that apparent R-wrapper defect; the Mobius score formulas are unchanged.

## Plotting

`dependogram` and `CrossCorrelogram` depend on `ggplot2` and are omitted.
The returned lag arrays and statistics can be written to a file or plotted
by a separate application.
