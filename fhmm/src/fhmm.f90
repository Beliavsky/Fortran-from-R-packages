! SPDX-License-Identifier: GPL-3.0-only
! Public facade for the numerical translation of fHMM 1.4.3.
module fhmm
   use fhmm_kinds
   use fhmm_types
   use fhmm_math, only: stationary_distribution, normal_cdf, normal_quantile
   use fhmm_distributions
   use fhmm_parameters
   use fhmm_algorithms
   use fhmm_hierarchical
   use fhmm_estimation
   use fhmm_diagnostics
   use fhmm_calendar
   implicit none
   public
end module fhmm
