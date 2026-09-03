! Copyright (C) 2000 Robert Gray
! Modern Fortran translation maintained for Fortran-from-R-packages.
! SPDX-License-Identifier: GPL-2.0-or-later
module cmprsk_summary
   use r_kinds, only : dp
   use r_distributions, only : r_pnorm, r_qnorm
   use cmprsk_status, only : cmprsk_success, cmprsk_invalid_argument
   use cmprsk_crr, only : crr_result
   implicit none
   private

   type, public :: crr_summary_result
      real(dp), allocatable :: coefficient(:)
      real(dp), allocatable :: relative_risk(:)
      real(dp), allocatable :: standard_error(:)
      real(dp), allocatable :: z_statistic(:)
      real(dp), allocatable :: p_value(:)
      real(dp), allocatable :: inverse_relative_risk(:)
      real(dp), allocatable :: confidence_lower(:)
      real(dp), allocatable :: confidence_upper(:)
      real(dp) :: likelihood_ratio = 0.0_dp
      integer :: degrees_freedom = 0
      real(dp) :: confidence_level = 0.95_dp
   end type crr_summary_result

   public :: summarize_crr

contains

   pure subroutine summarize_crr(fit, summary, status, confidence_level)
      type(crr_result), intent(in) :: fit !! Fine-Gray fit whose numerical summary statistics are requested.
      type(crr_summary_result), intent(out) :: summary !! Coefficient tests, relative risks, confidence limits, and LR test.
      integer, intent(out) :: status !! `cmprsk_success` or `cmprsk_invalid_argument`.
      real(dp), intent(in), optional :: confidence_level !! Two-sided level in `(0,1)`; defaults to `0.95`.

      integer :: i
      integer :: p
      real(dp) :: alpha
      real(dp) :: level
      real(dp) :: lower_z
      real(dp) :: upper_z

      level = 0.95_dp
      if (present(confidence_level)) level = confidence_level
      p = size(fit%coefficients)
      if (level <= 0.0_dp .or. level >= 1.0_dp .or. size(fit%variance, 1) /= p .or. &
          size(fit%variance, 2) /= p) then
         status = cmprsk_invalid_argument
         return
      end if
      if (any([(fit%variance(i, i), i = 1, p)] < 0.0_dp)) then
         status = cmprsk_invalid_argument
         return
      end if

      allocate(summary%coefficient(p), summary%relative_risk(p), summary%standard_error(p))
      allocate(summary%z_statistic(p), summary%p_value(p), summary%inverse_relative_risk(p))
      allocate(summary%confidence_lower(p), summary%confidence_upper(p))
      summary%coefficient = fit%coefficients
      summary%relative_risk = exp(fit%coefficients)
      summary%inverse_relative_risk = exp(-fit%coefficients)
      do i = 1, p
         summary%standard_error(i) = sqrt(fit%variance(i, i))
         if (summary%standard_error(i) > 0.0_dp) then
            summary%z_statistic(i) = fit%coefficients(i)/summary%standard_error(i)
            summary%p_value(i) = 2.0_dp*r_pnorm(abs(summary%z_statistic(i)), lower_tail=.false.)
         else
            summary%z_statistic(i) = 0.0_dp
            summary%p_value(i) = 1.0_dp
         end if
      end do
      alpha = 0.5_dp*(1.0_dp - level)
      lower_z = r_qnorm(alpha)
      upper_z = r_qnorm(1.0_dp - alpha)
      summary%confidence_lower = exp(fit%coefficients + lower_z*summary%standard_error)
      summary%confidence_upper = exp(fit%coefficients + upper_z*summary%standard_error)
      summary%likelihood_ratio = -2.0_dp*(fit%loglik_null - fit%loglik)
      summary%degrees_freedom = p
      summary%confidence_level = level
      status = cmprsk_success
   end subroutine summarize_crr

end module cmprsk_summary
