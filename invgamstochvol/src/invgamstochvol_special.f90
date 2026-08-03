! SPDX-License-Identifier: MIT
module invgamstochvol_special
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use invgamstochvol_kinds, only : dp
   use invgamstochvol_status, only : invgam_success, invgam_invalid_argument, &
      invgam_nonfinite_input, invgam_numerical_failure
   implicit none
   private

   public :: ourgeo, log_rising_factorial, build_log_factorials
   public :: hypergeo_from_tables, log_sum_exp

contains

   function log_rising_factorial(x, order, status) result(value)
      real(dp), intent(in) :: x
      integer, intent(in) :: order
      integer, intent(out), optional :: status
      real(dp) :: value
      integer :: j

      if (present(status)) status = invgam_success
      if (.not. ieee_is_finite(x) .or. order < 0) then
         value = 0.0_dp
         if (present(status)) status = invgam_invalid_argument
         return
      end if
      if (order == 0) then
         value = 0.0_dp
         return
      end if
      if (x <= 0.0_dp) then
         value = 0.0_dp
         if (present(status)) status = invgam_invalid_argument
         return
      end if

      value = 0.0_dp
      do j = 0, order - 1
         if (x + real(j, dp) <= 0.0_dp) then
            value = 0.0_dp
            if (present(status)) status = invgam_invalid_argument
            return
         end if
         value = value + log(x + real(j, dp))
      end do
   end function log_rising_factorial

   function ourgeo(a1, a2, b1, zstar, niter, status) result(value)
      real(dp), intent(in) :: a1, a2, b1, zstar
      integer, intent(in), optional :: niter
      integer, intent(out), optional :: status
      real(dp) :: value
      real(dp) :: term
      integer :: s, n_terms

      if (present(status)) status = invgam_success
      n_terms = 500
      if (present(niter)) n_terms = niter

      if (n_terms < 1 .or. .not. ieee_is_finite(a1) .or. &
          .not. ieee_is_finite(a2) .or. .not. ieee_is_finite(b1) .or. &
          .not. ieee_is_finite(zstar)) then
         value = 0.0_dp
         if (present(status)) status = invgam_invalid_argument
         return
      end if
      if (b1 <= 0.0_dp .or. zstar < 0.0_dp .or. zstar >= 1.0_dp) then
         value = 0.0_dp
         if (present(status)) status = invgam_invalid_argument
         return
      end if

      value = 1.0_dp
      term = 1.0_dp
      do s = 1, n_terms - 1
         term = term * (a1 + real(s - 1, dp)) * (a2 + real(s - 1, dp))
         term = term * zstar / ((b1 + real(s - 1, dp)) * real(s, dp))
         value = value + term
         if (.not. ieee_is_finite(value)) then
            if (present(status)) status = invgam_numerical_failure
            return
         end if
      end do
   end function ourgeo

   subroutine build_log_factorials(nu, nit, niter, alogfac, alogfac2, alfac, status)
      real(dp), intent(in) :: nu
      integer, intent(in) :: nit, niter
      real(dp), allocatable, intent(out) :: alogfac(:, :), alogfac2(:), alfac(:)
      integer, intent(out) :: status
      integer :: h, hold, max_order

      status = invgam_success
      if (nu <= 0.0_dp .or. nit < 0 .or. niter < 1) then
         status = invgam_invalid_argument
         return
      end if

      max_order = max(nit, niter)
      allocate(alogfac(0:nit, 0:max_order))
      allocate(alogfac2(0:max_order), alfac(0:max_order))

      do h = 0, max_order
         do hold = 0, nit
            alogfac(hold, h) = log_rising_factorial(0.5_dp * (nu + 1.0_dp) + &
               real(hold, dp), h, status)
            if (status /= invgam_success) return
         end do
         alogfac2(h) = log_rising_factorial(0.5_dp * nu, h, status)
         if (status /= invgam_success) return
         alfac(h) = log_rising_factorial(1.0_dp, h, status)
         if (status /= invgam_success) return
      end do
   end subroutine build_log_factorials

   function hypergeo_from_tables(h, alogfac, alogfac2, alfac, zstar, niter, status) result(value)
      integer, intent(in) :: h
      real(dp), intent(in) :: alogfac(0:, 0:), alogfac2(0:), alfac(0:)
      real(dp), intent(in) :: zstar
      integer, intent(in) :: niter
      integer, intent(out), optional :: status
      real(dp) :: value
      real(dp) :: log_z, log_term
      integer :: s

      if (present(status)) status = invgam_success
      if (h < 0 .or. h > ubound(alogfac, 1) .or. niter < 1 .or. &
          zstar < 0.0_dp .or. zstar >= 1.0_dp) then
         value = 0.0_dp
         if (present(status)) status = invgam_invalid_argument
         return
      end if
      if (niter - 1 > ubound(alogfac, 2)) then
         value = 0.0_dp
         if (present(status)) status = invgam_invalid_argument
         return
      end if

      value = 1.0_dp
      if (zstar <= tiny(1.0_dp)) return
      log_z = log(zstar)
      do s = 1, niter - 1
         log_term = alogfac(h, s) + alogfac(0, s) - alogfac2(s) &
            + real(s, dp) * log_z - alfac(s)
         value = value + exp(log_term)
      end do
      if (.not. ieee_is_finite(value)) then
         if (present(status)) status = invgam_numerical_failure
      end if
   end function hypergeo_from_tables

   function log_sum_exp(values, status) result(value)
      real(dp), intent(in) :: values(:)
      integer, intent(out), optional :: status
      real(dp) :: value
      real(dp) :: vmax, total
      integer :: i

      if (present(status)) status = invgam_success
      if (size(values) == 0) then
         value = -huge(1.0_dp)
         if (present(status)) status = invgam_invalid_argument
         return
      end if

      vmax = maxval(values)
      if (.not. ieee_is_finite(vmax)) then
         value = vmax
         if (present(status)) status = invgam_nonfinite_input
         return
      end if
      total = 0.0_dp
      do i = 1, size(values)
         total = total + exp(values(i) - vmax)
      end do
      if (total <= 0.0_dp .or. .not. ieee_is_finite(total)) then
         value = -huge(1.0_dp)
         if (present(status)) status = invgam_numerical_failure
         return
      end if
      value = vmax + log(total)
   end function log_sum_exp

end module invgamstochvol_special
