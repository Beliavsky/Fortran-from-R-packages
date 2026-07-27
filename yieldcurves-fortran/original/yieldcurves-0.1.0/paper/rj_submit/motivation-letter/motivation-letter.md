---
output: pdf_document
fontsize: 12pt
---

\thispagestyle{empty}
\today

The Editor
The R Journal
\bigskip

Dear Editor,
\bigskip

Please consider the article *Yieldcurves: Yield Curve Fitting, Analysis, and Decomposition in R* for publication in the R Journal.

The `yieldcurves` package consolidates the standard toolkit of term-structure analysis into a single CRAN package: Nelson-Siegel and Svensson parametric fitting, cubic-spline interpolation, extraction of forward rates and discount factors, Macaulay and modified duration, convexity, Z-spread, key rate durations, principal-component decomposition, level-slope-curvature factors, slope and butterfly measures, and carry-and-roll-down analysis. Existing CRAN infrastructure for this workflow is thin: `YieldCurve` fits Nelson-Siegel and Svensson but computes no risk measures and no factor decomposition, and the more comprehensive `termstrc` was archived in 2018. `yieldcurves` fills that gap with a uniform `yc_` prefixed API, plain numeric inputs, no heavy dependencies beyond `cli`, and a curve-object-based interface that lets downstream functions dispatch on a single S3 class.

Readers of the R Journal working in fixed-income research, quantitative asset management, central banking, and the teaching of financial markets should find the package useful. Applied uses include estimating a live Treasury curve, computing duration and convexity for a bond portfolio, running a PCA-based risk decomposition of a monthly panel of curves, and decomposing expected holding-period return into carry and roll-down contributions.

The manuscript has not been published in a peer-reviewed journal, is not currently under review elsewhere, and all rights to submit rest with the sole author.

\bigskip
\bigskip

Regards,
\bigskip
\bigskip

Charles Coverdale
London, United Kingdom
charles.f.coverdale@gmail.com
