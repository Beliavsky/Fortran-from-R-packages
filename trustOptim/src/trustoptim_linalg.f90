! SPDX-License-Identifier: MPL-2.0
module trustoptim_linalg
   use trustoptim_kinds, only : dp
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private
   public :: vecnorm, cholesky_lower, cholesky_solve, scaled_norm

contains

   pure function vecnorm(x) result(r)
      real(dp), intent(in) :: x(:)
      real(dp) :: r, scale, ssq, ax
      integer :: i

      scale = 0.0_dp
      ssq = 1.0_dp
      do i = 1, size(x)
         if (abs(x(i)) > 0.0_dp) then
            ax = abs(x(i))
            if (scale < ax) then
               ssq = 1.0_dp + ssq * (scale / ax)**2
               scale = ax
            else
               ssq = ssq + (ax / scale)**2
            end if
         end if
      end do
      if (scale <= 0.0_dp) then
         r = 0.0_dp
      else
         r = scale * sqrt(ssq)
      end if
   end function vecnorm

   subroutine cholesky_lower(a, l, ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: l(:,:)
      logical, intent(out) :: ok
      real(dp) :: s
      integer :: n, i, j, k

      n = size(a,1)
      l = 0.0_dp
      ok = .false.
      do j = 1, n
         s = a(j,j)
         do k = 1, j - 1
            s = s - l(j,k) * l(j,k)
         end do
         if (.not. (s > 0.0_dp) .or. .not. ieee_finite_local(s)) return
         l(j,j) = sqrt(s)
         do i = j + 1, n
            s = a(i,j)
            do k = 1, j - 1
               s = s - l(i,k) * l(j,k)
            end do
            l(i,j) = s / l(j,j)
         end do
      end do
      ok = .true.
   end subroutine cholesky_lower

   subroutine cholesky_solve(l, b, x)
      real(dp), intent(in) :: l(:,:), b(:)
      real(dp), intent(out) :: x(:)
      real(dp), allocatable :: y(:)
      integer :: n, i, j

      n = size(b)
      allocate(y(n))
      do i = 1, n
         y(i) = b(i)
         do j = 1, i - 1
            y(i) = y(i) - l(i,j) * y(j)
         end do
         y(i) = y(i) / l(i,i)
      end do
      do i = n, 1, -1
         x(i) = y(i)
         do j = i + 1, n
            x(i) = x(i) - l(j,i) * x(j)
         end do
         x(i) = x(i) / l(i,i)
      end do
   end subroutine cholesky_solve

   pure function scaled_norm(l, x) result(r)
      real(dp), intent(in) :: l(:,:), x(:)
      real(dp) :: r
      real(dp) :: y(size(x))
      integer :: i, j, n

      n = size(x)
      y = 0.0_dp
      do i = 1, n
         do j = i, n
            y(i) = y(i) + l(j,i) * x(j)
         end do
      end do
      r = vecnorm(y)
   end function scaled_norm

   pure logical function ieee_finite_local(x)
      real(dp), intent(in) :: x
      ieee_finite_local = ieee_is_finite(x)
   end function ieee_finite_local

end module trustoptim_linalg
