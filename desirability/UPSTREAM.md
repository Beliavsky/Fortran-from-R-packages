# Upstream provenance

This project is a source translation of the computational portions of the R package
`desirability`, version 2.1, authored by Max Kuhn.

Upstream package metadata is retained verbatim in `upstream/DESCRIPTION`. The upstream
package declares `License: GPL-2`; the complete GNU GPL version 2 text is in `COPYING`.

Translated source material:

- `R/desirability.R`: constructors and overall-combination logic
- `R/predict.R`: numerical prediction logic, interpolation, missing handling,
  tolerance handling, categorical handling, and geometric-mean aggregation

Intentionally not translated:

- `R/plot.R`: plotting/graphics code, per the translation scope
- `R/print.R`: R-specific presentation/S3 printing rather than computational logic
