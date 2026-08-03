! SPDX-License-Identifier: GPL-2.0-or-later
module ragtop_math
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use ragtop_kinds, only : dp, pi, sqrt_two
   use ragtop_constants, only : ragtop_ok, ragtop_invalid_argument, ragtop_singular_system
   implicit none
   private
   public :: normal_cdf, normal_pdf, linear_interp, cubic_spline_interp
   public :: solve_tridiagonal, sort_unique_real, bracketed_root

   abstract interface
      function scalar_function(x) result(y)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_function
   end interface

contains

   elemental real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp * erfc(-x / sqrt_two)
   end function normal_cdf

   elemental real(dp) function normal_pdf(x) result(p)
      real(dp), intent(in) :: x
      p = exp(-0.5_dp*x*x) / sqrt(2.0_dp*pi)
   end function normal_pdf

   pure real(dp) function linear_interp(x, y, xout) result(v)
      real(dp), intent(in) :: x(:), y(:), xout
      integer :: lo, hi, mid, n
      real(dp) :: w
      n = size(x)
      if (n == 0) then
         v = 0.0_dp
      else if (n == 1) then
         v = y(1)
      else if (xout <= x(1)) then
         w = (xout-x(1))/(x(2)-x(1))
         v = y(1) + w*(y(2)-y(1))
      else if (xout >= x(n)) then
         w = (xout-x(n-1))/(x(n)-x(n-1))
         v = y(n-1) + w*(y(n)-y(n-1))
      else
         lo = 1
         hi = n
         do while (hi-lo > 1)
            mid = (lo+hi)/2
            if (x(mid) <= xout) then
               lo = mid
            else
               hi = mid
            end if
         end do
         w = (xout-x(lo))/(x(hi)-x(lo))
         v = y(lo) + w*(y(hi)-y(lo))
      end if
   end function linear_interp

   function cubic_spline_interp(x, y, xout) result(v)
      real(dp), intent(in) :: x(:), y(:), xout
      real(dp) :: v
      real(dp), allocatable :: y2(:), u(:)
      real(dp) :: sig, p, h, a, b
      integer :: n, i, k_lo, k_hi
      n = size(x)
      if (n < 3) then
         v = linear_interp(x, y, xout)
         return
      end if
      allocate(y2(n), u(n-1))
      y2(1) = 0.0_dp
      u(1) = 0.0_dp
      do i = 2, n-1
         sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
         p = sig*y2(i-1) + 2.0_dp
         y2(i) = (sig-1.0_dp)/p
         u(i) = (6.0_dp*((y(i+1)-y(i))/(x(i+1)-x(i)) - &
                    (y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1)) - sig*u(i-1))/p
      end do
      y2(n) = 0.0_dp
      do i = n-1, 1, -1
         y2(i) = y2(i)*y2(i+1) + u(i)
      end do
      if (xout <= x(1) .or. xout >= x(n)) then
         v = linear_interp(x, y, xout)
         return
      end if
      k_lo = 1
      k_hi = n
      do while (k_hi-k_lo > 1)
         i = (k_hi+k_lo)/2
         if (x(i) > xout) then
            k_hi = i
         else
            k_lo = i
         end if
      end do
      h = x(k_hi)-x(k_lo)
      a = (x(k_hi)-xout)/h
      b = (xout-x(k_lo))/h
      v = a*y(k_lo) + b*y(k_hi) + &
          ((a*a*a-a)*y2(k_lo)+(b*b*b-b)*y2(k_hi))*h*h/6.0_dp
   end function cubic_spline_interp

   subroutine solve_tridiagonal(sub, diag, super, rhs, x, status)
      real(dp), intent(in) :: sub(:), diag(:), super(:), rhs(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: c(:), d(:)
      real(dp) :: denom
      integer :: n, i, st
      n = size(diag)
      st = ragtop_ok
      if (size(rhs) /= n .or. size(x) /= n .or. size(sub) /= max(0,n-1) .or. &
          size(super) /= max(0,n-1) .or. n < 1) then
         st = ragtop_invalid_argument
         x = 0.0_dp
         if (present(status)) status = st
         return
      end if
      allocate(c(max(1,n-1)), d(n))
      denom = diag(1)
      if (abs(denom) <= epsilon(denom)) then
         st = ragtop_singular_system
         x = 0.0_dp
         if (present(status)) status = st
         return
      end if
      if (n > 1) c(1) = super(1)/denom
      d(1) = rhs(1)/denom
      do i = 2, n
         denom = diag(i)-sub(i-1)*c(max(1,i-1))
         if (abs(denom) <= epsilon(denom)) then
            st = ragtop_singular_system
            x = 0.0_dp
            if (present(status)) status = st
            return
         end if
         if (i < n) c(i) = super(i)/denom
         d(i) = (rhs(i)-sub(i-1)*d(i-1))/denom
      end do
      x(n) = d(n)
      do i = n-1, 1, -1
         x(i) = d(i)-c(i)*x(i+1)
      end do
      if (present(status)) status = st
   end subroutine solve_tridiagonal

   subroutine sort_unique_real(x, y, tolerance)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable, intent(out) :: y(:)
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: work(:), tmp(:)
      real(dp) :: tol, key
      integer :: i, j, n, m
      n = size(x)
      tol = 100.0_dp*epsilon(1.0_dp)
      if (present(tolerance)) tol = max(0.0_dp,tolerance)
      allocate(work(n))
      work = x
      do i = 2, n
         key = work(i)
         j = i-1
         do while (j >= 1)
            if (work(j) <= key) exit
            work(j+1) = work(j)
            j = j-1
         end do
         work(j+1) = key
      end do
      if (n == 0) then
         allocate(y(0))
         return
      end if
      allocate(tmp(n))
      m = 1
      tmp(1) = work(1)
      do i = 2, n
         if (abs(work(i)-tmp(m)) > tol*max(1.0_dp,abs(work(i)),abs(tmp(m)))) then
            m = m+1
            tmp(m) = work(i)
         end if
      end do
      allocate(y(m))
      y = tmp(1:m)
   end subroutine sort_unique_real

   subroutine bracketed_root(f, target, lower, upper, root, status, tol, max_iter)
      procedure(scalar_function) :: f
      real(dp), intent(in) :: target, lower, upper
      real(dp), intent(out) :: root
      integer, intent(out) :: status
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: max_iter
      real(dp) :: a, b, c, fa, fb, fc, xtol
      integer :: iter, nmax
      xtol = 1.0e-8_dp
      if (present(tol)) xtol = max(tol,epsilon(1.0_dp))
      nmax = 100
      if (present(max_iter)) nmax = max(1,max_iter)
      a = lower
      b = upper
      fa = f(a)-target
      fb = f(b)-target
      if (ieee_is_nan(fa) .or. ieee_is_nan(fb) .or. fa*fb > 0.0_dp) then
         root = 0.5_dp*(a+b)
         status = ragtop_invalid_argument
         return
      end if
      do iter = 1, nmax
         c = 0.5_dp*(a+b)
         fc = f(c)-target
         if (abs(fc) <= xtol*max(1.0_dp,abs(target)) .or. abs(b-a) <= xtol*max(1.0_dp,abs(c))) then
            root = c
            status = ragtop_ok
            return
         end if
         if (fa*fc <= 0.0_dp) then
            b = c
            fb = fc
         else
            a = c
            fa = fc
         end if
      end do
      root = 0.5_dp*(a+b)
      status = 3
   end subroutine bracketed_root

end module ragtop_math
