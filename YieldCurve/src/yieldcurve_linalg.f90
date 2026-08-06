! SPDX-License-Identifier: GPL-2.0-or-later
module yieldcurve_linalg
   use yieldcurve_kinds, only : dp
   use yieldcurve_status, only : yc_success, yc_dimension_error, yc_rank_deficient
   implicit none
   private

   public :: least_squares

contains

   subroutine least_squares(a, b, x, residual, stat)
      real(dp), intent(in) :: a(:, :), b(:)
      real(dp), intent(out) :: x(:), residual(:)
      integer, intent(out) :: stat

      real(dp), allocatable :: q(:, :), r(:, :), v(:), y(:)
      real(dp) :: correction, scale, tolerance
      integer :: i, j, m, n

      m = size(a, 1)
      n = size(a, 2)
      stat = yc_success
      x = 0.0_dp
      residual = 0.0_dp

      if (size(b) /= m .or. size(x) /= n .or. size(residual) /= m .or. m < n) then
         stat = yc_dimension_error
         return
      end if

      allocate(q(m, n), r(n, n), v(m), y(n))
      q = 0.0_dp
      r = 0.0_dp
      scale = max(1.0_dp, maxval(abs(a)))
      tolerance = 100.0_dp * epsilon(1.0_dp) * real(max(m, n), dp) * scale

      do j = 1, n
         v = a(:, j)
         do i = 1, j - 1
            r(i, j) = dot_product(q(:, i), v)
            v = v - r(i, j) * q(:, i)
         end do

         ! A second pass substantially improves modified Gram-Schmidt when
         ! the two Svensson curvature columns are close to collinear.
         do i = 1, j - 1
            correction = dot_product(q(:, i), v)
            r(i, j) = r(i, j) + correction
            v = v - correction * q(:, i)
         end do

         r(j, j) = norm2(v)
         if (r(j, j) <= tolerance) then
            stat = yc_rank_deficient
            return
         end if
         q(:, j) = v / r(j, j)
      end do

      y = matmul(transpose(q), b)
      do i = n, 1, -1
         x(i) = y(i)
         if (i < n) x(i) = x(i) - dot_product(r(i, i + 1:n), x(i + 1:n))
         x(i) = x(i) / r(i, i)
      end do

      residual = b - matmul(a, x)
   end subroutine least_squares

end module yieldcurve_linalg
