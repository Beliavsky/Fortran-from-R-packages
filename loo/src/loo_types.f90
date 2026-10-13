! SPDX-License-Identifier: GPL-3.0-or-later
! Numerical translation of loo 2.10.1. See NOTICE.md.
module loo_types
   use r_kinds, only: dp
   use matrixstats, only: col_vars
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   implicit none
   private
   public :: dp, importance_sampling_result, loo_result, loo_nan, summarize_pointwise
   public :: log_likelihood_callback
   type :: importance_sampling_result
      real(dp), allocatable :: log_weights(:, :)
      real(dp), allocatable :: normalized_log_weights(:, :)
      real(dp), allocatable :: norm_const_log(:)
      real(dp), allocatable :: pareto_k(:)
      real(dp), allocatable :: n_eff(:)
      real(dp), allocatable :: r_eff(:)
      integer, allocatable :: tail_len(:)
      integer :: status = 1
      character(len = 4) :: method = ''
   end type importance_sampling_result
   type :: loo_result
      ! LOO columns: elpd, mcse_elpd, p_loo, looic, influence k.
      ! WAIC columns: elpd, p_waic, waic. ELPD columns: elpd, ic.
      real(dp), allocatable :: pointwise(:, :)
      ! Estimate and standard-error columns. LOO summarizes columns 1,3,4.
      real(dp), allocatable :: estimates(:, :)
      type(importance_sampling_result) :: importance
      integer :: n_draws = 0
      integer :: status = 1
      character(len = 5) :: method = ''
   end type loo_result
   abstract interface
      subroutine log_likelihood_callback(observation, values)
         import dp
         integer, intent(in) :: observation !! One-based observation index requested by the evaluator.
         real(dp), intent(out) :: values(:) !! One finite log likelihood per posterior draw, in fixed draw order.
      end subroutine log_likelihood_callback
   end interface
contains
   pure function loo_nan() result(x)
      real(dp) :: x
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function loo_nan

   pure function summarize_pointwise(x) result(summary)
      real(dp), intent(in) :: x(:, :) !! Rows are observations and columns are additive predictive quantities.
      real(dp) :: summary(size(x, 2), 2)
      summary(:, 1) = sum(x, dim = 1)
      summary(:, 2) = sqrt(real(size(x, 1), dp) * col_vars(x))
   end function summarize_pointwise
end module loo_types
