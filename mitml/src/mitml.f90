! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
! Public numerical facade for the modern Fortran translation of R package mitml 0.4-5.
module mitml
   use r_kinds, only : dp
   use mitml_cluster, only : cluster_means, cluster_means_matrix
   use mitml_convergence, only : gelman_rubin, moving_average, reduced_acf, sd_proportion
   use mitml_likelihood, only : gaussian_lm_loglik, gaussian_lmm_loglik
   use mitml_pool, only : pool_confint, pool_estimates
   use mitml_r2, only : intraclass_correlation, multilevel_r2
   use mitml_tests, only : d1_test, d2_test, d3_test, d4_test, test_linear_constraints, test_transformed_constraints
   use mitml_types, only : MITML_ERR_ARGUMENT, MITML_ERR_DIMENSION, MITML_ERR_LINALG, MITML_ERR_NUMERIC, MITML_OK, &
      mi_test_result, multilevel_r2_result, pooled_estimates
   implicit none
   private

   public :: dp
   public :: MITML_OK
   public :: MITML_ERR_DIMENSION
   public :: MITML_ERR_ARGUMENT
   public :: MITML_ERR_LINALG
   public :: MITML_ERR_NUMERIC
   public :: pooled_estimates
   public :: mi_test_result
   public :: multilevel_r2_result
   public :: pool_estimates
   public :: pool_confint
   public :: d1_test
   public :: d2_test
   public :: d3_test
   public :: d4_test
   public :: test_linear_constraints
   public :: test_transformed_constraints
   public :: cluster_means
   public :: cluster_means_matrix
   public :: multilevel_r2
   public :: intraclass_correlation
   public :: gelman_rubin
   public :: sd_proportion
   public :: reduced_acf
   public :: moving_average
   public :: gaussian_lm_loglik
   public :: gaussian_lmm_loglik

end module mitml
