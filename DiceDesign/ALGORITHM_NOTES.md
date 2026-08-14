# Algorithm notes

## Translation approach

DiceDesign 1.10 contains mostly R numerical code plus two C source files. The
computational algorithms were translated to native Fortran rather than wrapped
through R or C. Serialized R data used by the algorithms were decoded and
embedded as Fortran constants.

## Preserved computational data

The release contains the complete NOLH and NOLHDR libraries shipped in
`NOLHdesigns.rda` and `NOLHDRdesigns.rda`, including the 257-run designs through
29 dimensions. The Greenwood critical-value table from `sysdata.rda` is also
embedded.

## Deliberate corrections

Several apparent upstream implementation defects were corrected rather than
copied literally:

1. `dmaxDesign`: correlations involving a proposed point are explicitly reset to
   zero when the point is outside the variogram range. The R code could leave a
   stale correlation from the previous point location.
2. `unscaleDesign(..., uniformize=TRUE)`: the upstream branch refers to an
   uninitialized `x`; the Fortran routine performs the intended empirical
   quantile inverse using the supplied initial design.
3. ESE inner loops: `inner_iterations=N` performs exactly N proposals. The R
   `while(count <= inner_it)` structure performs one extra proposal.
4. ESE column cycling starts with the first design column rather than implicitly
   skipping it on the first proposal.
5. `meshRatio`: the extrema use nearest-neighbor distances for all design points;
   the R loop omits the final point.
6. RSS bounds are validated independently for lower and upper vectors rather
   than using the upstream conjunction that can let one malformed vector pass.

## LHS objective updates

The original package contains algebraic incremental-update helpers for L2/C2 and
phi-p criteria. The public Fortran optimizers preserve the same proposal and
acceptance algorithms but recompute the requested objective after each proposed
swap. This is computationally more expensive for large designs but substantially
simpler and numerically safer, especially for centered/boundary LHS values.

## RSS projection CDFs

`sumof2uniforms_cdf` and `sumof3uniforms_cdf` reproduce the arithmetic and the
`1e-12` degeneracy threshold in `src/C_rss.c` exactly. No defensive [0,1] clamp
is applied because the upstream C can leave that interval in extremely
ill-conditioned near-threshold cases; exact package compatibility was preferred.

## Random-number generation

R's global RNG calls are replaced by an explicit Park-Miller state object. This
makes Fortran calls reproducible and reentrant but does not reproduce R's random
stream for the same integer seed.

## Omitted interface/plotting code

Graphics in `mstCriteria`, `rss2d`, and `rss3d` are omitted. R list classes,
row/column names, formula/data-frame handling, and `xDRDN` symbolic name
formatting are interface concerns and are not reproduced. The corresponding
numeric computations are retained.
