# This file is part of PIQP-R.
#
# Copyright (c) 2023 piqp-r authors
#
# This source code is licensed under the BSD 2-Clause License found in the
# LICENSE file in the root directory of this source tree.

#' R Interface to PIQP Solver
#'
#' @description PIQP is an Proximal Interior Point Quadratic
#'   Programming solver, which can solve dense and sparse quadratic
#'   programs described in described in Schwan, Jiang, Kuhn, and Jones
#'   (2023) (<https://arxiv.org/abs/2304.00290>). Combining an
#'   infeasible interior point method with the proximal method of
#'   multipliers, the algorithm can handle ill-conditioned convex QP
#'   problems without the need for linear independence of the
#'   constraints. The solver is written in header only 'C++ 14'
#'   leveraging the Eigen library for vectorized linear algebra. For
#'   small dense problems, vectorized instructions and cache locality
#'   can be exploited more efficiently. Allocation free problem
#'   updates and re-solves are also provided.
#'
#' @name piqp-package
#' @useDynLib piqp
#' @importFrom Rcpp evalCpp
#' @author Balasubramanian Narasimhan, Roland Schwan (C), Yuning Jiang, Daniel Kuhn, Colin N. Jones
#' @keywords internal
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  S7::methods_register()
}
