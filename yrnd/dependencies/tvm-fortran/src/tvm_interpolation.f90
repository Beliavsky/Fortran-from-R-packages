! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Juan Manuel Truppia
! Modern Fortran translation of the tvm package.
module tvm_interpolation
   use tvm_kinds, only : dp
   implicit none
   private

   type, public :: pchip_interpolator
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: slope(:)
   contains
      procedure :: evaluate => pchip_evaluate_scalar
      procedure :: evaluate_many => pchip_evaluate_many
   end type pchip_interpolator

   public :: build_pchip

contains

   subroutine build_pchip(interp, x, y, status)
      type(pchip_interpolator), intent(out) :: interp
      real(dp), intent(in) :: x(:), y(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: h(:), delta(:)
      real(dp) :: w1, w2
      integer :: i, n, istat

      istat = 0
      n = size(x)
      if (n /= size(y) .or. n < 2) then
         istat = 1
         if (present(status)) status = istat
         return
      end if
      if (any(x(2:n) <= x(1:n - 1))) then
         istat = 2
         if (present(status)) status = istat
         return
      end if

      allocate(interp%x(n), interp%y(n), interp%slope(n))
      interp%x = x
      interp%y = y
      allocate(h(n - 1), delta(n - 1))
      h = x(2:n) - x(1:n - 1)
      delta = (y(2:n) - y(1:n - 1)) / h

      if (n == 2) then
         interp%slope = delta(1)
         if (present(status)) status = istat
         return
      end if

      interp%slope(1) = endpoint_slope(h(1), h(2), delta(1), delta(2))
      do i = 2, n - 1
         if (delta(i - 1) * delta(i) <= 0.0_dp) then
            interp%slope(i) = 0.0_dp
         else
            w1 = 2.0_dp * h(i) + h(i - 1)
            w2 = h(i) + 2.0_dp * h(i - 1)
            interp%slope(i) = (w1 + w2) / (w1 / delta(i - 1) + w2 / delta(i))
         end if
      end do
      interp%slope(n) = endpoint_slope(h(n - 1), h(n - 2), delta(n - 1), delta(n - 2))
      if (present(status)) status = istat
   end subroutine build_pchip

   pure real(dp) function endpoint_slope(h1, h2, delta1, delta2) result(slope)
      real(dp), intent(in) :: h1, h2, delta1, delta2

      slope = ((2.0_dp * h1 + h2) * delta1 - h1 * delta2) / (h1 + h2)
      if (slope * delta1 <= 0.0_dp) then
         slope = 0.0_dp
      else if (delta1 * delta2 < 0.0_dp .and. abs(slope) > 3.0_dp * abs(delta1)) then
         slope = 3.0_dp * delta1
      end if
   end function endpoint_slope

   pure real(dp) function pchip_evaluate_scalar(self, xq) result(yq)
      class(pchip_interpolator), intent(in) :: self
      real(dp), intent(in) :: xq
      real(dp) :: h, t, h00, h10, h01, h11
      integer :: i, n

      n = size(self%x)
      if (xq <= self%x(1)) then
         yq = self%y(1) + self%slope(1) * (xq - self%x(1))
         return
      end if
      if (xq >= self%x(n)) then
         yq = self%y(n) + self%slope(n) * (xq - self%x(n))
         return
      end if

      i = locate_interval(self%x, xq)
      h = self%x(i + 1) - self%x(i)
      t = (xq - self%x(i)) / h
      h00 = (1.0_dp + 2.0_dp * t) * (1.0_dp - t) ** 2
      h10 = t * (1.0_dp - t) ** 2
      h01 = t ** 2 * (3.0_dp - 2.0_dp * t)
      h11 = t ** 2 * (t - 1.0_dp)
      yq = h00 * self%y(i) + h10 * h * self%slope(i) + &
         h01 * self%y(i + 1) + h11 * h * self%slope(i + 1)
   end function pchip_evaluate_scalar

   pure function pchip_evaluate_many(self, xq) result(yq)
      class(pchip_interpolator), intent(in) :: self
      real(dp), intent(in) :: xq(:)
      real(dp) :: yq(size(xq))
      integer :: i

      do i = 1, size(xq)
         yq(i) = self%evaluate(xq(i))
      end do
   end function pchip_evaluate_many

   pure integer function locate_interval(x, xq) result(index)
      real(dp), intent(in) :: x(:), xq
      integer :: left, right, middle

      left = 1
      right = size(x)
      do while (right - left > 1)
         middle = (left + right) / 2
         if (xq < x(middle)) then
            right = middle
         else
            left = middle
         end if
      end do
      index = left
   end function locate_interval

end module tvm_interpolation
