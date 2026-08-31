! SPDX-License-Identifier: GPL-2.0-only
! Computational translation of CRAN mitools 2.4 by Thomas Lumley.
! Fortran translation and modifications: 2026-08-30.
module mitools_combine
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   use r_descriptive, only : r_covariance
   use r_distributions, only : r_qt
   use mitools_types, only : mi_result, mitools_insufficient_imputations
   use mitools_types, only : mitools_invalid_probability, mitools_invalid_shape, mitools_success
   implicit none
   private

   public :: mi_combine, mi_confidence_intervals, mi_standard_errors, mi_summary

   interface mi_combine
      module procedure mi_combine_matrix
      module procedure mi_combine_scalar
   end interface mi_combine

contains

   pure subroutine mi_combine_matrix(estimates, variances, result, status, df_complete)
      real(dp), intent(in) :: estimates(:, :) !! Parameter estimates, shape (n_parameter, n_imputation).
      real(dp), intent(in) :: variances(:, :, :) !! Within-imputation covariance matrices, shape (p, p, m).
      type(mi_result), intent(out) :: result !! Combined estimates, covariance, degrees of freedom, and missing information.
      integer, intent(out) :: status !! Zero on success; nonzero for incompatible shapes or fewer than two imputations.
      real(dp), intent(in), optional :: df_complete !! Complete-data degrees of freedom; omit for the upstream infinite default.
      real(dp), allocatable :: between(:, :)
      real(dp), allocatable :: within(:, :)
      real(dp) :: between_diag
      real(dp) :: complete_df
      real(dp) :: df_old
      real(dp) :: df_observed
      real(dp) :: relative_increase
      real(dp) :: within_diag
      integer :: i
      integer :: j
      integer :: m
      integer :: p

      status = mitools_success
      result%nimp = 0
      p = size(estimates, 1)
      m = size(estimates, 2)

      if (size(variances, 1) /= p .or. size(variances, 2) /= p .or. size(variances, 3) /= m) then
         status = mitools_invalid_shape
         return
      end if
      if (m < 2) then
         status = mitools_insufficient_imputations
         return
      end if

      complete_df = ieee_value(0.0_dp, ieee_positive_inf)
      if (present(df_complete)) complete_df = df_complete
      if (ieee_is_nan(complete_df) .or. complete_df <= 0.0_dp) then
         status = mitools_invalid_shape
         return
      end if

      allocate(result%coefficients(p), result%variance(p, p), result%df(p), result%missinfo(p))
      allocate(within(p, p), between(p, p))

      result%nimp = m
      result%coefficients = sum(estimates, dim=2) / real(m, dp)
      within = sum(variances, dim=3) / real(m, dp)

      do j = 1, p
         do i = 1, p
            between(i, j) = r_covariance(estimates(i, :), estimates(j, :))
         end do
      end do

      result%variance = within + between * (real(m + 1, dp) / real(m, dp))

      do i = 1, p
         within_diag = within(i, i)
         between_diag = between(i, i)
         relative_increase = relative_variance_increase(within_diag, between_diag, m)
         df_old = rubin_df(relative_increase, m)

         if (ieee_is_finite(complete_df)) then
            df_observed = ((complete_df + 1.0_dp) / (complete_df + 3.0_dp)) * complete_df * &
               within_diag / (within_diag + between_diag)
            result%df(i) = harmonic_df(df_old, df_observed)
         else
            result%df(i) = df_old
         end if

         result%missinfo(i) = (relative_increase + 2.0_dp / (result%df(i) + 3.0_dp)) / &
            (relative_increase + 1.0_dp)
      end do
   end subroutine mi_combine_matrix

   pure subroutine mi_combine_scalar(estimates, variances, result, status, df_complete)
      real(dp), intent(in) :: estimates(:) !! Scalar estimates from each imputation, length n_imputation.
      real(dp), intent(in) :: variances(:) !! Scalar within-imputation variances corresponding to the estimates.
      type(mi_result), intent(out) :: result !! Combined one-parameter multiple-imputation result.
      integer, intent(out) :: status !! Zero on success; nonzero for incompatible lengths or too few imputations.
      real(dp), intent(in), optional :: df_complete !! Complete-data degrees of freedom; omit for the upstream infinite default.
      real(dp), allocatable :: estimate_matrix(:, :)
      real(dp), allocatable :: variance_cube(:, :, :)
      integer :: m

      if (size(estimates) /= size(variances)) then
         status = mitools_invalid_shape
         result%nimp = 0
         return
      end if

      m = size(estimates)
      allocate(estimate_matrix(1, m), variance_cube(1, 1, m))
      estimate_matrix(1, :) = estimates
      variance_cube(1, 1, :) = variances
      if (present(df_complete)) then
         call mi_combine_matrix(estimate_matrix, variance_cube, result, status, df_complete)
      else
         call mi_combine_matrix(estimate_matrix, variance_cube, result, status)
      end if
   end subroutine mi_combine_scalar

   pure subroutine mi_standard_errors(result, standard_errors, status)
      type(mi_result), intent(in) :: result !! Combined multiple-imputation result whose covariance diagonal is summarized.
      real(dp), allocatable, intent(out) :: standard_errors(:) !! Square roots of the combined covariance diagonal.
      integer, intent(out) :: status !! Zero on success; nonzero when the result covariance has incompatible shape.
      integer :: i
      integer :: p

      status = mitools_success
      if (.not. allocated(result%coefficients) .or. .not. allocated(result%variance)) then
         status = mitools_invalid_shape
         return
      end if
      p = size(result%coefficients)
      if (size(result%variance, 1) /= p .or. size(result%variance, 2) /= p) then
         status = mitools_invalid_shape
         return
      end if

      allocate(standard_errors(p))
      do i = 1, p
         if (result%variance(i, i) < 0.0_dp) then
            standard_errors(i) = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            standard_errors(i) = sqrt(result%variance(i, i))
         end if
      end do
   end subroutine mi_standard_errors

   pure subroutine mi_confidence_intervals(result, alpha, lower, upper, status, log_effect)
      type(mi_result), intent(in) :: result !! Combined result supplying estimates, covariance, and Rubin degrees of freedom.
      real(dp), intent(in) :: alpha !! Two-sided error probability in the open interval (0, 1), such as 0.05.
      real(dp), allocatable, intent(out) :: lower(:) !! Lower confidence limits, optionally exponentiated for log effects.
      real(dp), allocatable, intent(out) :: upper(:) !! Upper confidence limits, optionally exponentiated for log effects.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid alpha or malformed result arrays.
      logical, intent(in), optional :: log_effect !! If true, exponentiate the confidence limits as in summary.MIresult.
      real(dp), allocatable :: standard_errors(:)
      real(dp) :: critical
      integer :: i
      integer :: p
      logical :: exponentiate

      status = mitools_success
      if (alpha <= 0.0_dp .or. alpha >= 1.0_dp .or. ieee_is_nan(alpha)) then
         status = mitools_invalid_probability
         return
      end if
      if (.not. allocated(result%coefficients) .or. .not. allocated(result%df)) then
         status = mitools_invalid_shape
         return
      end if
      p = size(result%coefficients)
      if (size(result%df) /= p) then
         status = mitools_invalid_shape
         return
      end if

      call mi_standard_errors(result, standard_errors, status)
      if (status /= mitools_success) return
      allocate(lower(p), upper(p))
      exponentiate = .false.
      if (present(log_effect)) exponentiate = log_effect

      do i = 1, p
         critical = r_qt(alpha / 2.0_dp, result%df(i), lower_tail=.false.)
         lower(i) = result%coefficients(i) - critical * standard_errors(i)
         upper(i) = result%coefficients(i) + critical * standard_errors(i)
         if (exponentiate) then
            lower(i) = exp(lower(i))
            upper(i) = exp(upper(i))
         end if
      end do
   end subroutine mi_confidence_intervals

   pure subroutine mi_summary(result, alpha, estimates, standard_errors, lower, upper, missinfo, status, log_effect)
      type(mi_result), intent(in) :: result !! Combined result to summarize in the style of summary.MIresult.
      real(dp), intent(in) :: alpha !! Two-sided error probability in the open interval (0, 1).
      real(dp), allocatable, intent(out) :: estimates(:) !! Point estimates, exponentiated when log_effect is true.
      real(dp), allocatable, intent(out) :: standard_errors(:) !! Standard errors, delta-method transformed for log effects.
      real(dp), allocatable, intent(out) :: lower(:) !! Lower confidence limits on the requested effect scale.
      real(dp), allocatable, intent(out) :: upper(:) !! Upper confidence limits on the requested effect scale.
      real(dp), allocatable, intent(out) :: missinfo(:) !! Fraction of missing information for each parameter.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid alpha or malformed result arrays.
      logical, intent(in), optional :: log_effect !! If true, report exponentiated effects and delta-method standard errors.
      integer :: p
      logical :: exponentiate

      call mi_standard_errors(result, standard_errors, status)
      if (status /= mitools_success) return
      call mi_confidence_intervals(result, alpha, lower, upper, status, log_effect)
      if (status /= mitools_success) return
      if (.not. allocated(result%missinfo)) then
         status = mitools_invalid_shape
         return
      end if

      p = size(result%coefficients)
      if (size(result%missinfo) /= p) then
         status = mitools_invalid_shape
         return
      end if
      allocate(estimates(p), missinfo(p))
      estimates = result%coefficients
      missinfo = result%missinfo

      exponentiate = .false.
      if (present(log_effect)) exponentiate = log_effect
      if (exponentiate) then
         estimates = exp(estimates)
         standard_errors = standard_errors * estimates
      end if
   end subroutine mi_summary

   pure elemental real(dp) function relative_variance_increase(within_variance, between_variance, m) result(value)
      real(dp), intent(in) :: within_variance !! Mean within-imputation variance for one parameter.
      real(dp), intent(in) :: between_variance !! Sample variance of that parameter across imputations.
      integer, intent(in) :: m !! Number of imputations, required to be at least two by the caller.

      if (ieee_is_nan(within_variance) .or. ieee_is_nan(between_variance)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (within_variance == 0.0_dp .and. between_variance == 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (within_variance == 0.0_dp .and. between_variance > 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_positive_inf)
      else
         value = (1.0_dp + 1.0_dp / real(m, dp)) * between_variance / within_variance
      end if
   end function relative_variance_increase

   pure elemental real(dp) function rubin_df(relative_increase, m) result(value)
      real(dp), intent(in) :: relative_increase !! Relative increase in variance from nonresponse for one parameter.
      integer, intent(in) :: m !! Number of imputations used in Rubin's degrees-of-freedom formula.

      if (ieee_is_nan(relative_increase)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (relative_increase == 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_positive_inf)
      else
         value = real(m - 1, dp) * (1.0_dp + 1.0_dp / relative_increase)**2
      end if
   end function rubin_df

   pure elemental real(dp) function harmonic_df(first_df, second_df) result(value)
      real(dp), intent(in) :: first_df !! First positive degrees-of-freedom contribution, typically Rubin's value.
      real(dp), intent(in) :: second_df !! Second positive contribution from the complete-data adjustment.

      if (ieee_is_nan(first_df) .or. ieee_is_nan(second_df)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (.not. ieee_is_finite(first_df)) then
         value = second_df
      else if (.not. ieee_is_finite(second_df)) then
         value = first_df
      else if (first_df == 0.0_dp .or. second_df == 0.0_dp) then
         value = 0.0_dp
      else
         value = 1.0_dp / (1.0_dp / first_df + 1.0_dp / second_df)
      end if
   end function harmonic_df

end module mitools_combine
