module countdm_math
   use countdm_kinds, only: dp
   implicit none
   private
   public :: lambert_w0, logistic, logit, logsumexp2, invert_matrix, numerical_hessian

contains

   pure real(dp) function logistic(x) result(p)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         p = 1.0_dp / (1.0_dp + exp(-x))
      else
         p = exp(x) / (1.0_dp + exp(x))
      end if
   end function logistic

   pure real(dp) function logit(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: q
      q = min(max(p, 1.0e-12_dp), 1.0_dp - 1.0e-12_dp)
      x = log(q / (1.0_dp - q))
   end function logit

   pure real(dp) function logsumexp2(a, b) result(v)
      real(dp), intent(in) :: a, b
      real(dp) :: m
      m = max(a, b)
      if (m < -0.5_dp * huge(1.0_dp)) then
         v = m
      else
         v = m + log(exp(a - m) + exp(b - m))
      end if
   end function logsumexp2

   pure real(dp) function lambert_w0(x) result(w)
      real(dp), intent(in) :: x
      real(dp) :: ew, f, den, dw
      integer :: i
      if (x < 0.0_dp) then
         w = ieee_nan()
         return
      else if (x == 0.0_dp) then
         w = 0.0_dp
         return
      end if
      if (x < 1.0_dp) then
         w = x
      else
         w = log(x)
         if (x > 3.0_dp) w = w - log(w)
      end if
      do i = 1, 50
         ew = exp(w)
         f = w * ew - x
         den = ew * (w + 1.0_dp) - (w + 2.0_dp) * f / (2.0_dp * w + 2.0_dp)
         dw = f / den
         w = w - dw
         if (abs(dw) <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(w))) exit
      end do
   contains
      pure real(dp) function ieee_nan() result(v)
         use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
         v = ieee_value(0.0_dp, ieee_quiet_nan)
      end function ieee_nan
   end function lambert_w0

   subroutine invert_matrix(a, ainv, ok)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: ainv(:, :)
      logical, intent(out) :: ok
      real(dp), allocatable :: aug(:, :), tmp(:)
      real(dp) :: pivot, factor
      integer :: n, i, j, k, p
      n = size(a, 1)
      ok = .false.
      if (size(a, 2) /= n) return
      allocate(aug(n, 2*n), ainv(n, n), tmp(2*n))
      aug = 0.0_dp
      aug(:, 1:n) = a
      do i = 1, n
         aug(i, n+i) = 1.0_dp
      end do
      do k = 1, n
         p = k
         do i = k + 1, n
            if (abs(aug(i, k)) > abs(aug(p, k))) p = i
         end do
         if (abs(aug(p, k)) <= 100.0_dp * epsilon(1.0_dp)) return
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
            aug(i, :) = aug(i, :) - factor * aug(k, :)
         end do
      end do
      ainv = aug(:, n+1:2*n)
      ok = .true.
   end subroutine invert_matrix

   subroutine numerical_hessian(fun, x, hess)
      interface
         function fun(z) result(f)
            import dp
            real(dp), intent(in) :: z(:)
            real(dp) :: f
         end function fun
      end interface
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: hess(:, :)
      real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:)
      real(dp) :: hi, hj, f0
      integer :: n, i, j
      n = size(x)
      allocate(xp(n), xm(n), xpp(n), xpm(n), xmp(n), xmm(n))
      f0 = fun(x)
      hess = 0.0_dp
      do i = 1, n
         hi = 2.0e-4_dp * max(1.0_dp, abs(x(i)))
         xp = x; xm = x
         xp(i) = xp(i) + hi
         xm(i) = xm(i) - hi
         hess(i, i) = (fun(xp) - 2.0_dp * f0 + fun(xm)) / (hi * hi)
         do j = i + 1, n
            hj = 2.0e-4_dp * max(1.0_dp, abs(x(j)))
            xpp = x; xpm = x; xmp = x; xmm = x
            xpp(i) = xpp(i) + hi; xpp(j) = xpp(j) + hj
            xpm(i) = xpm(i) + hi; xpm(j) = xpm(j) - hj
            xmp(i) = xmp(i) - hi; xmp(j) = xmp(j) + hj
            xmm(i) = xmm(i) - hi; xmm(j) = xmm(j) - hj
            hess(i, j) = (fun(xpp) - fun(xpm) - fun(xmp) + fun(xmm)) / (4.0_dp * hi * hj)
            hess(j, i) = hess(i, j)
         end do
      end do
   end subroutine numerical_hessian

end module countdm_math
