#' spantest: Mean–Variance Spanning Tests in R
#'
#' `spantest` implements classical and modern methods to test whether a set of
#' benchmark assets spans the mean–variance frontier of a larger asset universe.
#' The package includes Huberman–Kandel (1987), Gibbons–Ross–Shanken (1989),
#' Kempf–Memmel (2006), Pesaran–Yamagata (2024), Kan–Zhou (2012),
#' Britten-Jones (1999), and recent robust/high-dimensional tests such as
#' Gungor–Luger (2016) and Ardia–Sessinou (2025).
#'
#' @section Main functions:
#' - Alpha spanning: [span_grs()], [span_py()], [span_bj()], [span_f1()], [span_gl_a()], [span_as()]
#' - Variance spanning: [span_km()], [span_f2()], [span_as()]
#' - Joint mean–variance: [span_hk()], [span_gl_ad()], [span_as()]
#'
#' @section Documentation:
#' - Quick start: see the package README.
#' - Formal reference: the R Journal article (see `citation("spantest")` for the bib entry).
#' - Additional details and examples: function help pages and the README.
#'
#' @seealso
#' Related functions are grouped via `@family` tags for easy navigation.
#' Visit the package website for articles and examples.
#'
#' @references
#' \insertRef{ArdiaSessinou2025}{spantest} \cr
#' \insertRef{BrittenJones1999}{spantest} \cr
#' \insertRef{GRS1989}{spantest} \cr
#' \insertRef{GungorLuger2016}{spantest} \cr
#' \insertRef{HubermanKandel1987}{spantest} \cr
#' \insertRef{KanZhou2012}{spantest} \cr
#' \insertRef{KempfMemmel2006}{spantest} \cr
#' \insertRef{PesaranYamagata2024}{spantest} \cr
#'
#' @author David Ardia, Benjamin Seguin
#'
#' @name spantest
#' @keywords internal
#' @useDynLib spantest, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rdpack reprompt
"_PACKAGE"
