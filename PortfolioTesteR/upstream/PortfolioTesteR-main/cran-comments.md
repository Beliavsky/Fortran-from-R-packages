## Test environments
* Local: Windows 11, R 4.2.1 — `rcmdcheck --as-cran` => 0 errors | 0 warnings | 0 notes

## R CMD check results
0 errors | 0 warnings | 0 notes

## Release summary (0.1.3)
Maintenance/docs update; no breaking API changes.

### What changed
* Vignettes are built with **rmarkdown::html_vignette** and installed under **inst/doc**.
* Replaced non-ASCII symbols in Rd and vignettes to avoid LaTeX errors when building the PDF manual on Windows.
* Ensured no Quarto dependency is required for vignettes.

## CRAN policy notes
* Any examples that are long, networked, or heavy are wrapped in `\donttest{}` / guarded with `@examplesIf()`.
* Tests that use external resources call `testthat::skip_on_cran()`.
* Package does not write outside temp dirs; examples and vignettes use bundled data only.
* Declared data.table NSE globals to avoid NOTES.

## Reverse dependencies
None.
