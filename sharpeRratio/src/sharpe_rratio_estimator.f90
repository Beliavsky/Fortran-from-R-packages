! SPDX-License-Identifier: GPL-3.0-only
! Derived from sharpeRratio 1.4.3 by Damien Challet.
module sharpe_rratio_estimator
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ghyp_kinds, only : dp
   use ghyp_fitting, only : fit_result, fit_ghyp_uv
   use sharpe_rratio_records, only : r0_result, compute_r0bar
   use sharpe_rratio_statistics, only : test_n, theta_snr
   implicit none
   private

   real(dp), parameter, public :: gaussian_nu = 1.0e13_dp

   type, public :: snr_result
      real(dp) :: snr = 0.0_dp
      real(dp) :: ci_lower = 0.0_dp
      real(dp) :: ci_upper = 0.0_dp
      real(dp) :: nu = gaussian_nu
      real(dp) :: r0bar = 0.0_dp
      real(dp) :: normality_statistic = 0.0_dp
      integer :: n = 0
      integer :: num_permutations = 0
      logical :: nu_estimated = .false.
      logical :: gaussian_selected = .false.
      logical :: fit_converged = .false.
      logical :: ok = .false.
      character(len=200) :: message = ''
   end type snr_result

   public :: estimate_snr, estimateSNR, estimate_tail_exponent

   interface estimateSNR
      module procedure estimate_snr
   end interface estimateSNR

contains

   function estimate_tail_exponent(x, normality_statistic, max_fit_iterations, &
      fit_converged, gaussian_selected) result(nu)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out), optional :: normality_statistic
      integer, intent(in), optional :: max_fit_iterations
      logical, intent(out), optional :: fit_converged, gaussian_selected
      real(dp) :: nu, statistic
      type(fit_result) :: fit
      integer :: maxit
      logical :: gaussian

      statistic = test_n(x)
      if (present(normality_statistic)) normality_statistic = statistic
      gaussian = ieee_is_finite(statistic) .and. abs(statistic) < 3.0_dp
      if (gaussian) then
         nu = gaussian_nu
         if (present(fit_converged)) fit_converged = .true.
         if (present(gaussian_selected)) gaussian_selected = .true.
         return
      end if

      maxit = 1000
      if (present(max_fit_iterations)) maxit = max_fit_iterations
      fit = fit_ghyp_uv(x,'student',max_iter=maxit)
      if (fit%ok .and. fit%model%lambda < -1.0_dp) then
         nu = -2.0_dp*fit%model%lambda
      else
         nu = gaussian_nu
         gaussian = .true.
      end if
      if (present(fit_converged)) fit_converged = fit%converged .and. fit%ok
      if (present(gaussian_selected)) gaussian_selected = gaussian
   end function estimate_tail_exponent

   function estimate_snr(x, num_perm, nu, quantiles, seed, source_compatible, &
      max_fit_iterations) result(result)
      real(dp), intent(in) :: x(:)
      integer, intent(in), optional :: num_perm
      real(dp), intent(in), optional :: nu
      real(dp), intent(in), optional :: quantiles(2)
      integer(int64), intent(in), optional :: seed
      logical, intent(in), optional :: source_compatible
      integer, intent(in), optional :: max_fit_iterations
      type(snr_result) :: result
      real(dp), allocatable :: clean(:)
      real(dp) :: tail_exponent, probabilities(2)
      integer :: n, permutations
      logical :: fixed_nu, source_mode
      type(r0_result) :: records

      n = count(ieee_is_finite(x))
      if (n < 1) then
         result%message = 'at least one finite observation is required'
         return
      end if
      allocate(clean(n))
      clean = pack(x,ieee_is_finite(x))
      result%n = n

      permutations = min(100,max(1,ceiling(3.0_dp*log(real(n,dp)))))
      if (present(num_perm)) permutations = num_perm
      if (permutations < 1) then
         result%message = 'num_perm must be positive'
         return
      end if
      result%num_permutations = permutations

      fixed_nu = present(nu)
      if (fixed_nu) fixed_nu = nu >= 0.0_dp
      if (fixed_nu) then
         if (.not. ieee_is_finite(nu) .or. nu <= 0.0_dp) then
            result%message = 'nu must be positive and finite'
            return
         end if
         tail_exponent = nu
         result%normality_statistic = test_n(clean)
         result%fit_converged = .true.
      else
         tail_exponent = estimate_tail_exponent(clean,result%normality_statistic, &
            max_fit_iterations,result%fit_converged,result%gaussian_selected)
         result%nu_estimated = .true.
      end if
      result%nu = tail_exponent

      source_mode = .true.
      if (present(source_compatible)) source_mode = source_compatible
      if (source_mode) then
         probabilities = [0.025_dp,0.975_dp]
      else
         probabilities = [0.05_dp,0.95_dp]
         if (present(quantiles)) probabilities = quantiles
      end if
      if (any(probabilities < 0.0_dp) .or. any(probabilities > 1.0_dp)) then
         result%message = 'quantiles must lie in [0,1]'
         return
      end if

      records = compute_r0bar(clean,permutations,probabilities(1),probabilities(2), &
         seed,source_mode)
      if (.not. records%ok) then
         result%message = records%message
         return
      end if
      result%r0bar = records%mean
      if (abs(records%mean) < 1.0_dp) then
         result%snr = 0.0_dp
         result%ci_lower = 0.0_dp
         result%ci_upper = 0.0_dp
         result%ok = .true.
         return
      end if

      result%snr = theta_snr(records%mean/real(n,dp),n,tail_exponent,fixed_nu)
      result%ci_lower = theta_snr(records%q1/real(n,dp),n,tail_exponent,fixed_nu)
      result%ci_upper = theta_snr(records%q2/real(n,dp),n,tail_exponent,fixed_nu)
      if (result%ci_lower > result%ci_upper) then
         tail_exponent = result%ci_lower
         result%ci_lower = result%ci_upper
         result%ci_upper = tail_exponent
      end if
      result%ok = ieee_is_finite(result%snr) .and. &
         ieee_is_finite(result%ci_lower) .and. ieee_is_finite(result%ci_upper)
      if (.not. result%ok) result%message = 'non-finite estimate produced'
   end function estimate_snr

end module sharpe_rratio_estimator
