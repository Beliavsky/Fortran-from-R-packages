! SPDX-License-Identifier: GPL-2.0-only
module rsolnp_linalg
   use rsolnp_kinds, only : dp
   implicit none
   private

   public :: vnorm2, norm_inf, eye, outer_product
   public :: invert_matrix, symmetrize, projected_gradient

contains

   pure real(dp) function vnorm2(x) result(value)
      real(dp), intent(in) :: x(:)
      value = sqrt(max(0.0_dp, dot_product(x, x)))
   end function vnorm2

   pure real(dp) function norm_inf(x) result(value)
      real(dp), intent(in) :: x(:)
      if (size(x) == 0) then
         value = 0.0_dp
      else
         value = maxval(abs(x))
      end if
   end function norm_inf

   pure function eye(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n, n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end function eye

   pure function outer_product(x, y) result(a)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: a(size(x), size(y))
      integer :: i
      do i = 1, size(x)
         a(i, :) = x(i) * y
      end do
   end function outer_product

   pure subroutine symmetrize(a)
      real(dp), intent(inout) :: a(:, :)
      a = 0.5_dp * (a + transpose(a))
   end subroutine symmetrize

   subroutine invert_matrix(a, ainv, success)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(out) :: ainv(:, :)
      logical, intent(out) :: success

      real(dp), allocatable :: aug(:, :), tmp(:)
      real(dp) :: pivot, factor, scale
      integer :: n, i, k, p

      n = size(a, 1)
      success = .false.
      ainv = 0.0_dp
      if (size(a, 2) /= n .or. size(ainv, 1) /= n .or. size(ainv, 2) /= n) return
      allocate(aug(n, 2 * n), tmp(2 * n))
      aug(:, 1:n) = a
      aug(:, n + 1:2 * n) = eye(n)
      scale = max(1.0_dp, maxval(abs(a)))

      do k = 1, n
         p = k
         do i = k + 1, n
            if (abs(aug(i, k)) > abs(aug(p, k))) p = i
         end do
         if (abs(aug(p, k)) <= 100.0_dp * epsilon(1.0_dp) * scale) return
         if (p /= k) then
            tmp = aug(k, :)
            aug(k, :) = aug(p, :)
            aug(p, :) = tmp
         end if
         pivot = aug(k, k)
         aug(k, :) = aug(k, :) / pivot
         do i = 1, n
            if (i == k) cycle
            factor = aug(i, k)
            if (abs(factor) > tiny(1.0_dp)) aug(i, :) = aug(i, :) - factor * aug(k, :)
         end do
      end do
      ainv = aug(:, n + 1:2 * n)
      call symmetrize(ainv)
      success = .true.
   end subroutine invert_matrix

   pure subroutine projected_gradient(x, gradient, lower, upper, pg, active_tol)
      real(dp), intent(in) :: x(:), gradient(:), lower(:), upper(:)
      real(dp), intent(out) :: pg(:)
      real(dp), intent(in), optional :: active_tol
      real(dp) :: tol
      integer :: i

      tol = 100.0_dp * epsilon(1.0_dp)
      if (present(active_tol)) tol = active_tol
      pg = gradient
      do i = 1, size(x)
         if (x(i) <= lower(i) + tol * max(1.0_dp, abs(lower(i)))) then
            if (gradient(i) > 0.0_dp) pg(i) = 0.0_dp
         end if
         if (x(i) >= upper(i) - tol * max(1.0_dp, abs(upper(i)))) then
            if (gradient(i) < 0.0_dp) pg(i) = 0.0_dp
         end if
      end do
   end subroutine projected_gradient

end module rsolnp_linalg
