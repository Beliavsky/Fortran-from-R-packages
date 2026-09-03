module tmb_numerics
   use tmb_kinds, only: dp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
   implicit none
   private
   public :: kronecker_product, find_interval, sort_real, order_real, interpolate2d
   public :: cholesky_lower, mvnorm_nll, ar1_nll, ar1_mvn_nll, matrix_exponential, romberg_integrate
   public :: unstructured_corr, unstructured_corr_nll, n01_nll
   abstract interface
      pure function scalar_integrand(x) result(value)
         import dp
         real(dp), intent(in) :: x !! Scalar integration coordinate.
         real(dp) :: value
      end function scalar_integrand
   end interface
contains

   pure real(dp) function romberg_integrate(f, a, b, levels) result(ans)
      procedure(scalar_integrand) :: f !! Pure scalar integrand.
      real(dp), intent(in) :: a !! Lower integration limit.
      real(dp), intent(in) :: b !! Upper integration limit.
      integer, intent(in), optional :: levels !! Romberg levels; default 7 and must be at least 1.
      integer :: lev, i, j, n, m
      real(dp) :: h, subtotal
      real(dp), allocatable :: r(:, :)
      m = 7
      if (present(levels)) m = levels
      if (m < 1) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      allocate (r(m, m))
      r = 0.0_dp
      h = b - a
      r(1, 1) = 0.5_dp * h * (f(a) + f(b))
      do lev = 2, m
         n = 2**(lev - 2)
         subtotal = 0.0_dp
         do i = 1, n
            subtotal = subtotal + f(a + (real(i, dp) - 0.5_dp) * h)
         end do
         r(lev, 1) = 0.5_dp * r(lev - 1, 1) + 0.5_dp * h * subtotal
         do j = 2, lev
            r(lev, j) = r(lev, j - 1) + (r(lev, j - 1) - r(lev - 1, j - 1)) / (4.0_dp**(j - 1) - 1.0_dp)
         end do
         h = 0.5_dp * h
      end do
      ans = r(m, m)
   end function romberg_integrate

   pure function kronecker_product(a, b) result(c)
      real(dp), intent(in) :: a(:, :) !! Left matrix of shape (m,n).
      real(dp), intent(in) :: b(:, :) !! Right matrix of shape (p,q).
      real(dp) :: c(size(a, 1) * size(b, 1), size(a, 2) * size(b, 2))
      integer :: i, j, p, q
      p = size(b, 1)
      q = size(b, 2)
      do j = 1, size(a, 2)
         do i = 1, size(a, 1)
            c((i - 1) * p + 1:i * p, (j - 1) * q + 1:j * q) = a(i, j) * b
         end do
      end do
   end function kronecker_product

   pure integer function find_interval(x, breaks) result(idx)
      real(dp), intent(in) :: x !! Scalar query point.
      real(dp), intent(in) :: breaks(:) !! Nondecreasing breakpoints.
      integer :: lo, hi, mid
      if (size(breaks) == 0 .or. x < breaks(1)) then
         idx = 0
         return
      end if
      if (x >= breaks(size(breaks))) then
         idx = size(breaks)
         return
      end if
      lo = 1
      hi = size(breaks)
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (x < breaks(mid)) then
            hi = mid
         else
            lo = mid
         end if
      end do
      idx = lo
   end function find_interval

   pure function order_real(x) result(ord)
      real(dp), intent(in) :: x(:) !! Values to order in ascending order.
      integer :: ord(size(x))
      integer :: i, j, key
      ord = [(i, i=1, size(x))]
      do i = 2, size(x)
         key = ord(i)
         j = i - 1
         do while (j >= 1)
            if (x(ord(j)) <= x(key)) exit
            ord(j + 1) = ord(j)
            j = j - 1
         end do
         ord(j + 1) = key
      end do
   end function order_real

   pure function sort_real(x) result(y)
      real(dp), intent(in) :: x(:) !! Values to sort in ascending order.
      real(dp) :: y(size(x))
      integer :: ord(size(x))
      ord = order_real(x)
      y = x(ord)
   end function sort_real

   pure real(dp) function interpolate2d(data, x_range, y_range, x, y, radius) result(ans)
      real(dp), intent(in) :: data(:, :) !! Gridded data; NaNs are ignored.
      real(dp), intent(in) :: x_range(2) !! Coordinates of first and last matrix rows.
      real(dp), intent(in) :: y_range(2) !! Coordinates of first and last matrix columns.
      real(dp), intent(in) :: x !! First-coordinate query point.
      real(dp), intent(in) :: y !! Second-coordinate query point.
      real(dp), intent(in) :: radius !! Positive smoothing radius in lattice-index units.
      real(dp) :: hx, hy, xi, yi, dist, w, fw_sum, w_sum, pi
      integer :: i, j, imin, imax, jmin, jmax
      pi = acos(-1.0_dp)
      if (radius <= 0.0_dp .or. size(data, 1) < 2 .or. size(data, 2) < 2) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      hx = (x_range(2) - x_range(1)) / real(size(data, 1) - 1, dp)
      hy = (y_range(2) - y_range(1)) / real(size(data, 2) - 1, dp)
      xi = (x - x_range(1)) / hx
      yi = (y - y_range(1)) / hy
      if (xi < 0.0_dp .or. xi > real(size(data, 1) - 1, dp) .or. &
          yi < 0.0_dp .or. yi > real(size(data, 2) - 1, dp)) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      imin = max(0, ceiling(xi - radius))
      imax = min(size(data, 1) - 1, floor(xi + radius))
      jmin = max(0, ceiling(yi - radius))
      jmax = min(size(data, 2) - 1, floor(yi + radius))
      fw_sum = 0.0_dp
      w_sum = 0.0_dp
      do j = jmin, jmax
         do i = imin, imax
            dist = sqrt((real(i, dp) - xi)**2 + (real(j, dp) - yi)**2 + 1.0e-100_dp)
            if (dist <= radius .and. .not. ieee_is_nan(data(i + 1, j + 1))) then
               w = 0.5_dp * (1.0_dp + cos((1.0_dp - 0.5_dp * (1.0_dp + cos(pi * dist / radius))) * pi))
               fw_sum = fw_sum + data(i + 1, j + 1) * w
               w_sum = w_sum + w
            end if
         end do
      end do
      if (w_sum > 0.0_dp) then
         ans = fw_sum / w_sum
      else
         ans = ieee_value(ans, ieee_quiet_nan)
      end if
   end function interpolate2d

   pure subroutine cholesky_lower(a, l, info)
      real(dp), intent(in) :: a(:, :) !! Symmetric positive-definite input matrix.
      real(dp), intent(out) :: l(:, :) !! Lower-triangular Cholesky factor when info=0.
      integer, intent(out) :: info !! Zero on success; positive pivot index on failure.
      integer :: i, j, k, n
      real(dp) :: s
      n = size(a, 1)
      l = 0.0_dp
      info = 0
      if (size(a, 2) /= n .or. size(l, 1) /= n .or. size(l, 2) /= n) then
         info = -1
         return
      end if
      do j = 1, n
         s = a(j, j)
         do k = 1, j - 1
            s = s - l(j, k) * l(j, k)
         end do
         if (s <= 0.0_dp) then
            info = j
            return
         end if
         l(j, j) = sqrt(s)
         do i = j + 1, n
            s = a(i, j)
            do k = 1, j - 1
               s = s - l(i, k) * l(j, k)
            end do
            l(i, j) = s / l(j, j)
         end do
      end do
   end subroutine cholesky_lower

   pure real(dp) function mvnorm_nll(x, sigma) result(ans)
      real(dp), intent(in) :: x(:) !! Zero-mean Gaussian observation vector.
      real(dp), intent(in) :: sigma(:, :) !! Symmetric positive-definite covariance matrix.
      real(dp) :: l(size(x), size(x)), y(size(x)), s
      integer :: i, j, info, n
      n = size(x)
      if (size(sigma, 1) /= n .or. size(sigma, 2) /= n) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      call cholesky_lower(sigma, l, info)
      if (info /= 0) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      y = 0.0_dp
      do i = 1, n
         s = x(i)
         do j = 1, i - 1
            s = s - l(i, j) * y(j)
         end do
         y(i) = s / l(i, i)
      end do
      ans = 0.5_dp * real(n, dp) * log(2.0_dp * acos(-1.0_dp)) + sum(log([(l(i, i), i=1, n)])) + &
            0.5_dp * dot_product(y, y)
   end function mvnorm_nll

   pure elemental real(dp) function n01_nll(x) result(ans)
      real(dp), intent(in) :: x !! Standard-normal scalar observation.
      ans = 0.5_dp * x * x + 0.5_dp * log(2.0_dp * acos(-1.0_dp))
   end function n01_nll

   pure subroutine unstructured_corr(theta, sigma, info)
      real(dp), intent(in) :: theta(:) !! Strict-lower-triangle parameters filled row-wise.
      real(dp), intent(out) :: sigma(:, :) !! Correlation matrix implied by the unit-diagonal lower factor.
      integer, intent(out) :: info !! Zero on success; nonzero when theta length or sigma shape is invalid.
      real(dp), allocatable :: l(:, :), llt(:, :)
      real(dp) :: root
      integer :: i, j, k, m, n
      m = size(theta)
      root = sqrt(1.0_dp + 8.0_dp * real(m, dp))
      n = nint(0.5_dp * (1.0_dp + root))
      if (n * (n - 1) / 2 /= m .or. size(sigma, 1) /= n .or. size(sigma, 2) /= n) then
         sigma = ieee_value(0.0_dp, ieee_quiet_nan)
         info = 1
         return
      end if
      allocate (l(n, n), llt(n, n))
      l = 0.0_dp
      do i = 1, n
         l(i, i) = 1.0_dp
      end do
      k = 0
      do i = 1, n
         do j = 1, i - 1
            k = k + 1
            l(i, j) = theta(k)
         end do
      end do
      llt = matmul(l, transpose(l))
      do j = 1, n
         do i = 1, n
            sigma(i, j) = llt(i, j) / sqrt(llt(i, i) * llt(j, j))
         end do
      end do
      info = 0
   end subroutine unstructured_corr

   pure real(dp) function unstructured_corr_nll(x, theta) result(ans)
      real(dp), intent(in) :: x(:) !! Zero-mean observation vector.
      real(dp), intent(in) :: theta(:) !! Strict-lower-triangle correlation parameters.
      real(dp), allocatable :: sigma(:, :)
      integer :: info, n
      n = size(x)
      if (size(theta) /= n * (n - 1) / 2) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      allocate (sigma(n, n))
      call unstructured_corr(theta, sigma, info)
      if (info /= 0) then
         ans = ieee_value(ans, ieee_quiet_nan)
      else
         ans = mvnorm_nll(x, sigma)
      end if
   end function unstructured_corr_nll

   pure real(dp) function ar1_mvn_nll(x, phi, sigma) result(ans)
      real(dp), intent(in) :: x(:, :) !! Multivariate AR(1) states; columns are successive times.
      real(dp), intent(in) :: phi !! Scalar autoregressive coefficient strictly between -1 and 1.
      real(dp), intent(in) :: sigma(:, :) !! Marginal covariance matrix of each state vector.
      real(dp) :: innovation(size(x, 1)), scale
      integer :: i
      if (abs(phi) >= 1.0_dp .or. size(x, 2) == 0 .or. size(sigma, 1) /= size(x, 1) .or. &
          size(sigma, 2) /= size(x, 1)) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      scale = sqrt(1.0_dp - phi * phi)
      ans = mvnorm_nll(x(:, 1), sigma)
      do i = 2, size(x, 2)
         innovation = (x(:, i) - phi * x(:, i - 1)) / scale
         ans = ans + mvnorm_nll(innovation, sigma) + real(size(x, 1), dp) * log(scale)
      end do
   end function ar1_mvn_nll

   pure real(dp) function ar1_nll(x, phi) result(ans)
      real(dp), intent(in) :: x(:) !! Stationary AR(1) realization with unit marginal variance.
      real(dp), intent(in) :: phi !! Autoregressive coefficient strictly between -1 and 1.
      real(dp) :: innovation_sd, pi
      integer :: i
      pi = acos(-1.0_dp)
      if (abs(phi) >= 1.0_dp .or. size(x) == 0) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      innovation_sd = sqrt(1.0_dp - phi * phi)
      ans = 0.5_dp * log(2.0_dp * pi) + 0.5_dp * x(1) * x(1)
      do i = 2, size(x)
         ans = ans + log(innovation_sd) + 0.5_dp * log(2.0_dp * pi) + &
               0.5_dp * ((x(i) - phi * x(i - 1)) / innovation_sd)**2
      end do
   end function ar1_nll

   pure function matrix_exponential(a) result(ea)
      real(dp), intent(in) :: a(:, :) !! Square matrix whose exponential is requested.
      real(dp) :: ea(size(a, 1), size(a, 2))
      real(dp) :: b(size(a, 1), size(a, 2)), term(size(a, 1), size(a, 2))
      real(dp) :: norm1
      integer :: i, k, n, scale_pow
      n = size(a, 1)
      ea = 0.0_dp
      if (size(a, 2) /= n) then
         ea = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      norm1 = maxval(sum(abs(a), dim=1))
      if (norm1 > 0.5_dp) then
         scale_pow = max(0, ceiling(log(norm1 / 0.5_dp) / log(2.0_dp)))
      else
         scale_pow = 0
      end if
      b = a / 2.0_dp**scale_pow
      do i = 1, n
         ea(i, i) = 1.0_dp
      end do
      term = ea
      do k = 1, 40
         term = matmul(term, b) / real(k, dp)
         ea = ea + term
         if (maxval(abs(term)) < 1.0e-15_dp) exit
      end do
      do k = 1, scale_pow
         ea = matmul(ea, ea)
      end do
   end function matrix_exponential
end module tmb_numerics
