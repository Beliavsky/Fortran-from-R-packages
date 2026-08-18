! BiasedUrn-fortran
! Computational translation of the CRAN BiasedUrn package.
! Upstream copyright (c) 2002-2024 Agner Fog.
! License: GNU General Public License version 3.
module biasedurn_math
   use biasedurn_kinds, only : dp
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
   implicit none
   private

   real(dp), parameter, public :: log_zero = -huge(1.0_dp) / 4.0_dp
   real(dp), parameter :: pi = acos(-1.0_dp)

   public :: log_choose, log_add, clamp_precision
   public :: quiet_nan, positive_inf
   public :: biasedurn_seed, rand_uniform, sample_log_weights
   public :: wallenius_log_integral

contains

   pure function quiet_nan() result(x)
      real(dp) :: x
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

   pure function positive_inf() result(x)
      real(dp) :: x
      x = ieee_value(0.0_dp, ieee_positive_inf)
   end function positive_inf

   pure function log_choose(n, k) result(v)
      integer, intent(in) :: n, k
      real(dp) :: v
      if (k < 0 .or. k > n .or. n < 0) then
         v = log_zero
      else
         v = log_gamma(real(n + 1, dp)) - log_gamma(real(k + 1, dp)) &
            - log_gamma(real(n - k + 1, dp))
      end if
   end function log_choose

   pure function log_add(a, b) result(c)
      real(dp), intent(in) :: a, b
      real(dp) :: c, hi, lo
      if (a <= log_zero / 2.0_dp) then
         c = b
      else if (b <= log_zero / 2.0_dp) then
         c = a
      else
         hi = max(a, b)
         lo = min(a, b)
         c = hi + log1p_std(exp(lo - hi))
      end if
   end function log_add

   pure function clamp_precision(precision) result(p)
      real(dp), intent(in), optional :: precision
      real(dp) :: p
      if (present(precision)) then
         p = precision
      else
         p = 1.0e-7_dp
      end if
      if (.not. (p > 0.0_dp .and. p < 1.0_dp)) p = 1.0e-7_dp
      p = max(p, 1.0e-13_dp)
   end function clamp_precision

   subroutine biasedurn_seed(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)
      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(abs(seed) + 104729 * i + 8191 * i * i, huge(1) - 1)
         if (put(i) == 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine biasedurn_seed

   function rand_uniform() result(u)
      real(dp) :: u
      call random_number(u)
      if (u <= 0.0_dp) u = tiny(1.0_dp)
   end function rand_uniform

   integer function sample_log_weights(logw) result(index)
      real(dp), intent(in) :: logw(:)
      real(dp) :: m, s, u, c
      integer :: i

      m = maxval(logw)
      if (m <= log_zero / 2.0_dp) then
         index = 0
         return
      end if
      s = 0.0_dp
      do i = 1, size(logw)
         if (logw(i) > log_zero / 2.0_dp) s = s + exp(logw(i) - m)
      end do
      u = rand_uniform() * s
      c = 0.0_dp
      do i = 1, size(logw)
         if (logw(i) > log_zero / 2.0_dp) c = c + exp(logw(i) - m)
         if (u <= c) then
            index = i
            return
         end if
      end do
      index = size(logw)
   end function sample_log_weights

   function wallenius_log_integral(x, m, odds, precision) result(logi)
      ! Computes
      ! integral_0^1 product_i (1-t**(odds_i/d))**x_i dt,
      ! d = sum_i odds_i * (m_i-x_i), using t=exp(-u) and
      ! mode-scaled Gauss-Legendre integration.
      integer, intent(in) :: x(:), m(:)
      real(dp), intent(in) :: odds(:)
      real(dp), intent(in), optional :: precision
      real(dp) :: logi
      real(dp) :: d, p, mode, lmode, target, left, right, step
      real(dp) :: wmax
      real(dp), allocatable :: wn(:)
      real(dp) :: il, ir, a, b, fa, fb, mid
      integer :: i, iter

      if (size(x) /= size(m) .or. size(x) /= size(odds)) then
         logi = quiet_nan()
         return
      end if
      if (sum(x) == 0) then
         logi = 0.0_dp
         return
      end if

      wmax = maxval(odds)
      if (wmax <= 0.0_dp) then
         logi = quiet_nan()
         return
      end if
      allocate(wn(size(odds)))
      wn = odds / wmax
      d = 0.0_dp
      do i = 1, size(x)
         d = d + wn(i) * real(m(i) - x(i), dp)
      end do
      if (d <= 0.0_dp) then
         ! This occurs only at the deterministic endpoint n=sum(m)
         ! for valid positive-weight inputs. The defining integral tends to
         ! reciprocal multinomial coefficient there; callers handle it.
         logi = quiet_nan()
         return
      end if

      p = clamp_precision(precision)

      ! The derivative of log integrand is monotone decreasing:
      ! -1 + sum x_i*a_i/(exp(a_i*u)-1), a_i=odds_i/d.
      a = max(tiny(1.0_dp), 1.0e-14_dp)
      b = 1.0_dp
      do while (log_integrand_derivative(b, x, wn, d) > 0.0_dp)
         b = 2.0_dp * b
         if (b > 1.0e8_dp) exit
      end do
      do iter = 1, 100
         mid = 0.5_dp * (a + b)
         if (log_integrand_derivative(mid, x, wn, d) > 0.0_dp) then
            a = mid
         else
            b = mid
         end if
         if (abs(b - a) <= 1.0e-12_dp * max(1.0_dp, mid)) exit
      end do
      mode = 0.5_dp * (a + b)
      lmode = log_integrand(mode, x, wn, d)

      ! Truncate where the scaled integrand is safely below requested error.
      target = log(p) - 8.0_dp

      a = 0.0_dp
      b = mode
      do iter = 1, 100
         mid = 0.5_dp * (a + b)
         if (log_integrand(mid, x, wn, d) - lmode < target) then
            a = mid
         else
            b = mid
         end if
      end do
      left = b

      step = max(1.0_dp, 0.25_dp * max(1.0_dp, mode))
      right = mode + step
      do while (log_integrand(right, x, wn, d) - lmode > target)
         step = 1.7_dp * step
         right = mode + step
         if (right > mode + 1.0e6_dp) exit
      end do
      a = mode
      b = right
      do iter = 1, 100
         mid = 0.5_dp * (a + b)
         if (log_integrand(mid, x, wn, d) - lmode > target) then
            a = mid
         else
            b = mid
         end if
      end do
      right = b

      il = 0.0_dp
      if (mode > left) il = gl_integral_scaled(left, mode, x, wn, d, lmode)
      ir = gl_integral_scaled(mode, right, x, wn, d, lmode)
      if (il + ir <= 0.0_dp) then
         logi = log_zero
      else
         logi = lmode + log(il + ir)
      end if
   end function wallenius_log_integral

   pure function log_integrand(u, x, odds, d) result(v)
      real(dp), intent(in) :: u, odds(:), d
      integer, intent(in) :: x(:)
      real(dp) :: v, z, ai
      integer :: i
      if (u <= 0.0_dp) then
         if (sum(x) > 0) then
            v = log_zero
         else
            v = 0.0_dp
         end if
         return
      end if
      v = -u
      do i = 1, size(x)
         if (x(i) == 0) cycle
         ai = odds(i) / d
         if (ai <= 0.0_dp) then
            v = log_zero
            return
         end if
         z = ai * u
         if (z < 0.6931471805599453_dp) then
            v = v + real(x(i), dp) * log(-expm1_std(-z))
         else
            v = v + real(x(i), dp) * log1p_std(-exp(-z))
         end if
      end do
   end function log_integrand

   pure function log_integrand_derivative(u, x, odds, d) result(v)
      real(dp), intent(in) :: u, odds(:), d
      integer, intent(in) :: x(:)
      real(dp) :: v, z, ai, den
      integer :: i
      v = -1.0_dp
      do i = 1, size(x)
         if (x(i) == 0) cycle
         ai = odds(i) / d
         z = ai * u
         if (z > 700.0_dp) cycle
         den = expm1_std(z)
         if (den <= 0.0_dp) then
            v = huge(1.0_dp)
            return
         end if
         v = v + real(x(i), dp) * ai / den
      end do
   end function log_integrand_derivative

   function gl_integral_scaled(a, b, x, odds, d, lmode) result(s)
      real(dp), intent(in) :: a, b, odds(:), d, lmode
      integer, intent(in) :: x(:)
      real(dp) :: s
      real(dp) :: nodes(32), weights(32), xm, xr, u
      integer :: i
      call gauss_legendre_32(nodes, weights)
      xm = 0.5_dp * (a + b)
      xr = 0.5_dp * (b - a)
      s = 0.0_dp
      do i = 1, 32
         u = xm + xr * nodes(i)
         s = s + weights(i) * exp(log_integrand(u, x, odds, d) - lmode)
      end do
      s = xr * s
   end function gl_integral_scaled

   subroutine gauss_legendre_32(x, w)
      real(dp), intent(out) :: x(32), w(32)
      real(dp) :: z, z1, p1, p2, p3, pp
      integer :: i, j, m
      m = 16
      do i = 1, m
         z = cos(pi * (real(i, dp) - 0.25_dp) / 32.5_dp)
         do
            p1 = 1.0_dp
            p2 = 0.0_dp
            do j = 1, 32
               p3 = p2
               p2 = p1
               p1 = ((2.0_dp * real(j, dp) - 1.0_dp) * z * p2 &
                  - (real(j, dp) - 1.0_dp) * p3) / real(j, dp)
            end do
            pp = 32.0_dp * (z * p1 - p2) / (z * z - 1.0_dp)
            z1 = z
            z = z1 - p1 / pp
            if (abs(z - z1) <= 4.0_dp * epsilon(1.0_dp)) exit
         end do
         x(i) = -z
         x(33 - i) = z
         w(i) = 2.0_dp / ((1.0_dp - z * z) * pp * pp)
         w(33 - i) = w(i)
      end do
   end subroutine gauss_legendre_32

   pure function log1p_std(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: v, term, sum, old
      integer :: k
      if (abs(x) > 1.0e-4_dp) then
         v = log(1.0_dp + x)
         return
      end if
      sum = 0.0_dp
      term = x
      do k = 1, 200
         old = sum
         sum = sum + merge(1.0_dp, -1.0_dp, mod(k, 2) == 1) * term / real(k, dp)
         if (sum == old) exit
         term = term * x
      end do
      v = sum
   end function log1p_std

   pure function expm1_std(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: v, term, sum, old
      integer :: k
      if (abs(x) > 1.0e-5_dp) then
         v = exp(x) - 1.0_dp
         return
      end if
      sum = x
      term = x
      do k = 2, 100
         old = sum
         term = term * x / real(k, dp)
         sum = sum + term
         if (sum == old) exit
      end do
      v = sum
   end function expm1_std

end module biasedurn_math
