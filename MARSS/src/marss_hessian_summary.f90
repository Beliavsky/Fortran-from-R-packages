! SPDX-License-Identifier: GPL-2.0-only
module marss_hessian_summary_mod
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use marss_kinds, only : dp
   use marss_types, only : marss_model, marss_constraints, marss_hessian_result
   use marss_constraints_mod, only : marss_constraints_parameter_count
   use marss_constraints_mod, only : marss_fisher_i_linear, marss_vectorized_parameter_names
   use marss_harvey_info, only : marss_fisher_i_harvey_linear
   use r_linalg, only : inverse_matrix
   implicit none
   private
   public :: marss_hessian_summary_linear

contains

   subroutine marss_hessian_summary_linear(template, constraints, theta, result, method, rel_step)
      type(marss_model), intent(in) :: template !! Model containing data and fixed numerical blocks for information evaluation.
      type(marss_constraints), intent(in) :: constraints !! Affine constraints defining the estimated beta coordinates.
      real(dp), intent(in) :: theta(:) !! Fitted beta coordinates at which observed information is evaluated.
      type(marss_hessian_result), intent(out) :: result !! Parameter mean, observed information, covariance, names, and status.
      character(len=*), intent(in), optional :: method !! Information method: Harvey1989, fdHess, or optim; default is Harvey1989.
      real(dp), intent(in), optional :: rel_step !! Relative finite-difference step used by fdHess and optim modes.
      real(dp), allocatable :: work(:, :)
      logical, allocatable :: nan_diag(:)
      character(len=16) :: selected
      real(dp) :: nan_value
      integer :: i
      integer :: p
      integer :: stat

      result%info = 0
      result%covariance_info = 0
      result%covariance_available = .false.
      result%had_nan_information = .false.
      selected = 'harvey1989'
      if (present(method)) selected = canonical_method(method)
      result%method = selected

      p = marss_constraints_parameter_count(constraints)
      if (size(theta) /= p .or. p == 0) then
         result%info = 1
         return
      end if

      allocate(result%par_mean(p))
      result%par_mean = theta
      call marss_vectorized_parameter_names(constraints, result%parameter_names)

      select case (trim(selected))
      case ('harvey1989')
         call marss_fisher_i_harvey_linear(template, constraints, theta, result%hessian, stat)
      case ('fdhess', 'optim')
         if (present(rel_step)) then
            call marss_fisher_i_linear(template, constraints, theta, result%hessian, rel_step, stat)
         else
            call marss_fisher_i_linear(template, constraints, theta, result%hessian, info=stat)
         end if
      case default
         result%info = 2
         return
      end select
      if (stat /= 0) then
         result%info = 10 + stat
         return
      end if

      allocate(work(p, p), nan_diag(p))
      work = result%hessian
      nan_diag = .false.
      do i = 1, p
         nan_diag(i) = ieee_is_nan(work(i, i))
      end do
      result%had_nan_information = any(ieee_is_nan(work))
      do i = 1, p
         if (nan_diag(i)) work(i, i) = 1.0_dp
      end do
      where (ieee_is_nan(work)) work = 0.0_dp

      call inverse_matrix(work, result%par_sigma, stat)
      if (stat /= 0) then
         result%covariance_info = stat
         if (allocated(result%par_sigma)) deallocate(result%par_sigma)
         return
      end if

      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      do i = 1, p
         if (nan_diag(i)) result%par_sigma(i, i) = nan_value
      end do
      result%covariance_available = .true.
   end subroutine marss_hessian_summary_linear

   pure function canonical_method(text) result(value)
      character(len=*), intent(in) :: text !! User method spelling normalized without changing supported method semantics.
      character(len=16) :: value
      integer :: code
      integer :: i
      integer :: n

      value = ''
      n = min(len_trim(text), len(value))
      do i = 1, n
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            value(i:i) = achar(code + iachar('a') - iachar('A'))
         else
            value(i:i) = text(i:i)
         end if
      end do
      value = adjustl(value)
   end function canonical_method

end module marss_hessian_summary_mod
