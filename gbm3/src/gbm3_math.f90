! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of gbm3 3.0.3; see NOTICE.md and upstream/.
module gbm3_math
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use gbm3_kinds, only : dp
   implicit none
   private
   public :: sigmoid, softplus, log_add_exp, exp_diff, weighted_quantile, location_m_tdist
   public :: argsort_real, argsort_desc_real, shuffle_int, quiet_nan, safe_divide

contains

   pure real(dp) function quiet_nan() result(x)
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

   pure real(dp) function sigmoid(x) result(p)
      real(dp), intent(in) :: x
      real(dp) :: e
      if (x >= 0.0_dp) then
         e = exp(-x)
         p = 1.0_dp / (1.0_dp + e)
      else
         e = exp(x)
         p = e / (1.0_dp + e)
      end if
   end function sigmoid

   pure real(dp) function softplus(x) result(y)
      real(dp), intent(in) :: x
      if (x > 0.0_dp) then
         y = x + log(1.0_dp + exp(-x))
      else
         y = log(1.0_dp + exp(x))
      end if
   end function softplus

   pure real(dp) function log_add_exp(total, value) result(out)
      real(dp), intent(in) :: total, value
      if (total <= -huge(1.0_dp)) then
         out = value
      else if (value <= -huge(1.0_dp)) then
         out = total
      else if (value > total) then
         out = value + log(1.0_dp + exp(total - value))
      else
         out = total + log(1.0_dp + exp(value - total))
      end if
   end function log_add_exp

   pure real(dp) function exp_diff(log_x, log_y) result(out)
      real(dp), intent(in) :: log_x, log_y
      real(dp) :: m
      if (abs(log_x - log_y) <= tiny(1.0_dp)) then
         out = 0.0_dp
      else
         m = max(log_x, log_y)
         out = exp(m) * (exp(log_x - m) - exp(log_y - m))
      end if
   end function exp_diff

   pure real(dp) function safe_divide(num, den) result(out)
      real(dp), intent(in) :: num, den
      if (abs(den) <= tiny(1.0_dp)) then
         out = 0.0_dp
      else
         out = num / den
      end if
   end function safe_divide

   subroutine argsort_real(x, idx)
      real(dp), intent(in) :: x(:)
      integer, intent(out) :: idx(size(x))
      integer :: i, j, key
      real(dp) :: keyval
      do i = 1, size(x)
         idx(i) = i
      end do
      do i = 2, size(x)
         key = idx(i)
         keyval = x(key)
         j = i - 1
         do while (j >= 1)
            if (x(idx(j)) <= keyval) exit
            idx(j + 1) = idx(j)
            j = j - 1
         end do
         idx(j + 1) = key
      end do
   end subroutine argsort_real

   subroutine argsort_desc_real(x, idx)
      real(dp), intent(in) :: x(:)
      integer, intent(out) :: idx(size(x))
      integer :: i, j, key
      real(dp) :: keyval
      do i = 1, size(x)
         idx(i) = i
      end do
      do i = 2, size(x)
         key = idx(i)
         keyval = x(key)
         j = i - 1
         do while (j >= 1)
            if (x(idx(j)) >= keyval) exit
            idx(j + 1) = idx(j)
            j = j - 1
         end do
         idx(j + 1) = key
      end do
   end subroutine argsort_desc_real

   subroutine shuffle_int(x)
      integer, intent(inout) :: x(:)
      integer :: i, j, tmp
      real(dp) :: u
      do i = size(x), 2, -1
         call random_number(u)
         j = 1 + int(u * real(i, dp))
         if (j > i) j = i
         tmp = x(i)
         x(i) = x(j)
         x(j) = tmp
      end do
   end subroutine shuffle_int

   real(dp) function weighted_quantile(values, weights, alpha) result(q)
      real(dp), intent(in) :: values(:), weights(:), alpha
      integer, allocatable :: ord(:)
      real(dp) :: wsum, cum
      integer :: i, med_idx, next_nonzero

      if (size(values) /= size(weights)) error stop "weighted_quantile: shape mismatch"
      if (size(values) == 0) then
         q = 0.0_dp
         return
      else if (size(values) == 1) then
         q = values(1)
         return
      end if

      allocate(ord(size(values)))
      call argsort_real(values, ord)
      wsum = sum(weights)
      med_idx = 0
      cum = 0.0_dp
      do while (cum < alpha * wsum .and. med_idx < size(values))
         med_idx = med_idx + 1
         cum = cum + weights(ord(med_idx))
      end do
      if (med_idx == 0) med_idx = 1

      next_nonzero = size(values) + 1
      do i = size(values), med_idx + 1, -1
         if (weights(ord(i)) > 0.0_dp) next_nonzero = i
      end do

      if (next_nonzero == size(values) + 1 .or. cum > alpha * wsum) then
         q = values(ord(med_idx))
      else
         q = (1.0_dp - alpha) * values(ord(med_idx)) + alpha * values(ord(next_nonzero))
      end if
   end function weighted_quantile

   real(dp) function location_m_tdist(values, weights, nu, alpha) result(beta0)
      real(dp), intent(in) :: values(:), weights(:), nu, alpha
      real(dp), allocatable :: diff(:)
      real(dp) :: scale0, dt, dwt, sum_w_x, sum_w, beta, err
      real(dp), parameter :: meps = 1.0e-6_dp
      integer :: i, counter

      if (size(values) /= size(weights)) error stop "location_m_tdist: shape mismatch"
      if (size(values) == 0) then
         beta0 = 0.0_dp
         return
      end if

      beta0 = weighted_quantile(values, weights, alpha)
      allocate(diff(size(values)))
      diff = abs(values - beta0)
      scale0 = max(1.4826_dp * weighted_quantile(diff, weights, alpha), meps)

      do counter = 1, 50
         sum_w_x = 0.0_dp
         sum_w = 0.0_dp
         do i = 1, size(values)
            dt = max(abs(values(i) - beta0) / scale0, meps)
            dwt = weights(i) * (dt / (nu + dt * dt)) / dt
            sum_w_x = sum_w_x + dwt * values(i)
            sum_w = sum_w + dwt
         end do
         beta = beta0
         if (sum_w > 0.0_dp) beta = sum_w_x / sum_w
         err = abs(beta - beta0)
         if (err > meps .and. abs(beta0) > tiny(1.0_dp)) err = err / abs(beta0)
         beta0 = beta
         if (err < meps) exit
      end do
   end function location_m_tdist

end module gbm3_math
