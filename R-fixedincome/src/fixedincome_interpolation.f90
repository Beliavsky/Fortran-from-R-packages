! SPDX-License-Identifier: MIT
module fixedincome_interpolation
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use fixedincome_kinds, only : dp
   use fixedincome_types, only : interpolation_t, fit_result_t, &
      INTERP_FLAT_FORWARD, INTERPOLATION_LINEAR, INTERP_LOG_LINEAR, INTERP_NATURAL_SPLINE, &
      INTERP_HERMITE_SPLINE, INTERP_MONOTONE_SPLINE, INTERP_NELSON_SIEGEL, &
      INTERP_NELSON_SIEGEL_SVENSSON, FI_OK, FI_INVALID_ARGUMENT, FI_OUT_OF_RANGE, &
      FI_NO_CONVERGENCE
   implicit none
   private

   public :: interp_flatforward, interp_linear, interp_loglinear
   public :: interp_naturalspline, interp_hermitespline, interp_monotonespline
   public :: interp_nelsonsiegel, interp_nelsonsiegelsvensson
   public :: nelson_siegel, nelson_siegel_svensson, fit_interpolation_model
   public :: interpolate_points, parameters, interpolation_name

contains


   function parameters(model) result(values)
      type(interpolation_t), intent(in) :: model
      real(dp), allocatable :: values(:)
      select case (model%method)
      case (INTERP_NELSON_SIEGEL)
         values = model%parameters(:4)
      case (INTERP_NELSON_SIEGEL_SVENSSON)
         values = model%parameters(:6)
      case default
         allocate(values(0))
      end select
   end function parameters

   pure function interpolation_name(model) result(name)
      type(interpolation_t), intent(in) :: model
      character(len=32) :: name
      select case (model%method)
      case (INTERP_FLAT_FORWARD)
         name = 'flatforward'
      case (INTERPOLATION_LINEAR)
         name = 'linear'
      case (INTERP_LOG_LINEAR)
         name = 'loglinear'
      case (INTERP_NATURAL_SPLINE)
         name = 'naturalspline'
      case (INTERP_HERMITE_SPLINE)
         name = 'hermitespline'
      case (INTERP_MONOTONE_SPLINE)
         name = 'monotonespline'
      case (INTERP_NELSON_SIEGEL)
         name = 'nelsonsiegel'
      case (INTERP_NELSON_SIEGEL_SVENSSON)
         name = 'nelsonsiegelsvensson'
      case default
         name = 'none'
      end select
   end function interpolation_name

   pure function interp_flatforward(propagate) result(model)
      logical, intent(in), optional :: propagate
      type(interpolation_t) :: model
      model%method = INTERP_FLAT_FORWARD
      if (present(propagate)) model%propagate = propagate
   end function interp_flatforward

   pure function interp_linear(propagate) result(model)
      logical, intent(in), optional :: propagate
      type(interpolation_t) :: model
      model%method = INTERPOLATION_LINEAR
      if (present(propagate)) model%propagate = propagate
   end function interp_linear

   pure function interp_loglinear(propagate) result(model)
      logical, intent(in), optional :: propagate
      type(interpolation_t) :: model
      model%method = INTERP_LOG_LINEAR
      if (present(propagate)) model%propagate = propagate
   end function interp_loglinear

   pure function interp_naturalspline(propagate) result(model)
      logical, intent(in), optional :: propagate
      type(interpolation_t) :: model
      model%method = INTERP_NATURAL_SPLINE
      if (present(propagate)) model%propagate = propagate
   end function interp_naturalspline

   pure function interp_hermitespline(propagate) result(model)
      logical, intent(in), optional :: propagate
      type(interpolation_t) :: model
      model%method = INTERP_HERMITE_SPLINE
      if (present(propagate)) model%propagate = propagate
   end function interp_hermitespline

   pure function interp_monotonespline(propagate) result(model)
      logical, intent(in), optional :: propagate
      type(interpolation_t) :: model
      model%method = INTERP_MONOTONE_SPLINE
      if (present(propagate)) model%propagate = propagate
   end function interp_monotonespline

   pure function interp_nelsonsiegel(beta1, beta2, beta3, lambda1, propagate) result(model)
      real(dp), intent(in) :: beta1, beta2, beta3, lambda1
      logical, intent(in), optional :: propagate
      type(interpolation_t) :: model
      model%method = INTERP_NELSON_SIEGEL
      model%parameters(1:4) = [beta1, beta2, beta3, lambda1]
      if (present(propagate)) model%propagate = propagate
   end function interp_nelsonsiegel

   pure function interp_nelsonsiegelsvensson(beta1, beta2, beta3, beta4, lambda1, lambda2, &
                                              propagate) result(model)
      real(dp), intent(in) :: beta1, beta2, beta3, beta4, lambda1, lambda2
      logical, intent(in), optional :: propagate
      type(interpolation_t) :: model
      model%method = INTERP_NELSON_SIEGEL_SVENSSON
      model%parameters = [beta1, beta2, beta3, beta4, lambda1, lambda2]
      if (present(propagate)) model%propagate = propagate
   end function interp_nelsonsiegelsvensson

   elemental pure real(dp) function nelson_siegel(t, beta1, beta2, beta3, lambda1) result(rate)
      real(dp), intent(in) :: t, beta1, beta2, beta3, lambda1
      real(dp) :: loading1, loading2, z
      z = lambda1 * t
      if (abs(z) < sqrt(epsilon(1.0_dp))) then
         loading1 = 1.0_dp - 0.5_dp * z + z*z / 6.0_dp
      else
         loading1 = -expm1_safe(-z) / z
      end if
      loading2 = loading1 - exp(-z)
      rate = beta1 + beta2 * loading1 + beta3 * loading2
   end function nelson_siegel

   elemental pure real(dp) function nelson_siegel_svensson(t, beta1, beta2, beta3, beta4, &
                                                            lambda1, lambda2) result(rate)
      real(dp), intent(in) :: t, beta1, beta2, beta3, beta4, lambda1, lambda2
      real(dp) :: z2, loading
      z2 = lambda2 * t
      if (abs(z2) < sqrt(epsilon(1.0_dp))) then
         loading = (1.0_dp - 0.5_dp*z2 + z2*z2/6.0_dp) - exp(-z2)
      else
         loading = -expm1_safe(-z2) / z2 - exp(-z2)
      end if
      rate = nelson_siegel(t, beta1, beta2, beta3, lambda1) + beta4 * loading
   end function nelson_siegel_svensson

   function interpolate_points(method, x, y, query, parameters, status) result(values)
      integer, intent(in) :: method
      real(dp), intent(in) :: x(:), y(:), query(:)
      real(dp), intent(in), optional :: parameters(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: slopes(:), second(:)
      integer :: i, stat_i
      allocate(values(size(query)))
      if (method == INTERP_NELSON_SIEGEL .or. method == INTERP_NELSON_SIEGEL_SVENSSON) then
         if (.not. present(parameters)) then
            values = ieee_value(0.0_dp, ieee_quiet_nan)
            if (present(status)) status = FI_INVALID_ARGUMENT
            return
         end if
         select case (method)
         case (INTERP_NELSON_SIEGEL)
            if (size(parameters) < 4) then
               values = ieee_value(0.0_dp, ieee_quiet_nan)
               if (present(status)) status = FI_INVALID_ARGUMENT
               return
            end if
            values = nelson_siegel(query, parameters(1), parameters(2), parameters(3), parameters(4))
         case (INTERP_NELSON_SIEGEL_SVENSSON)
            if (size(parameters) < 6) then
               values = ieee_value(0.0_dp, ieee_quiet_nan)
               if (present(status)) status = FI_INVALID_ARGUMENT
               return
            end if
            values = nelson_siegel_svensson(query, parameters(1), parameters(2), parameters(3), &
                                             parameters(4), parameters(5), parameters(6))
         end select
         if (present(status)) status = FI_OK
         return
      end if

      if (size(x) /= size(y) .or. size(x) < 2 .or. any(x(2:) <= x(:size(x)-1))) then
         values = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end if

      select case (method)
      case (INTERPOLATION_LINEAR, INTERP_LOG_LINEAR, INTERP_FLAT_FORWARD)
         do i = 1, size(query)
            if (method == INTERP_LOG_LINEAR) then
               if (any(y <= 0.0_dp)) then
                  values = ieee_value(0.0_dp, ieee_quiet_nan)
                  if (present(status)) status = FI_INVALID_ARGUMENT
                  return
               end if
               values(i) = exp(linear_one(x, log(y), query(i), stat_i))
            else
               values(i) = linear_one(x, y, query(i), stat_i)
            end if
            if (stat_i /= FI_OK) then
               values(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
         end do
      case (INTERP_NATURAL_SPLINE)
         second = natural_second_derivatives(x, y, stat_i)
         if (stat_i /= FI_OK) then
            values = ieee_value(0.0_dp, ieee_quiet_nan)
            if (present(status)) status = stat_i
            return
         end if
         do i = 1, size(query)
            values(i) = natural_spline_one(x, y, second, query(i), stat_i)
            if (stat_i /= FI_OK) values(i) = ieee_value(0.0_dp, ieee_quiet_nan)
         end do
      case (INTERP_HERMITE_SPLINE, INTERP_MONOTONE_SPLINE)
         slopes = monotone_slopes(x, y)
         do i = 1, size(query)
            values(i) = hermite_one(x, y, slopes, query(i), stat_i)
            if (stat_i /= FI_OK) values(i) = ieee_value(0.0_dp, ieee_quiet_nan)
         end do
      case default
         values = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = FI_INVALID_ARGUMENT
         return
      end select
      if (present(status)) then
         if (any(query < x(1)) .or. any(query > x(size(x)))) then
            status = FI_OUT_OF_RANGE
         else
            status = FI_OK
         end if
      end if
   end function interpolate_points

   function fit_interpolation_model(initial_model, terms_years, rates, max_iterations, tolerance) result(fit)
      type(interpolation_t), intent(in) :: initial_model
      real(dp), intent(in) :: terms_years(:), rates(:)
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(fit_result_t) :: fit
      real(dp), allocatable :: simplex(:, :), values(:), centroid(:), xr(:), xe(:), xc(:)
      real(dp), allocatable :: lower_bound(:), upper_bound(:), step(:)
      real(dp) :: alpha, gamma, rho, sigma, tol, simplex_spread, fnew
      integer :: n, i, iter, maxit, best, worst, second_worst

      if (size(terms_years) /= size(rates) .or. size(rates) < 2) then
         fit%status = FI_INVALID_ARGUMENT
         return
      end if
      select case (initial_model%method)
      case (INTERP_NELSON_SIEGEL)
         n = 4
         allocate(lower_bound(n), upper_bound(n), step(n))
         lower_bound = [0.0_dp, -0.3_dp, -1.0_dp, 1.0e-6_dp]
         upper_bound = [0.3_dp, 0.3_dp, 1.0_dp, 5.0_dp]
         step = [0.01_dp, 0.02_dp, 0.03_dp, 0.1_dp]
      case (INTERP_NELSON_SIEGEL_SVENSSON)
         n = 6
         allocate(lower_bound(n), upper_bound(n), step(n))
         lower_bound = [0.0_dp, -0.3_dp, -1.0_dp, -1.0_dp, 1.0e-6_dp, 1.0e-6_dp]
         upper_bound = [0.3_dp, 0.3_dp, 1.0_dp, 1.0_dp, 5.0_dp, 5.0_dp]
         step = [0.01_dp, 0.02_dp, 0.03_dp, 0.03_dp, 0.1_dp, 0.1_dp]
      case default
         fit%status = FI_INVALID_ARGUMENT
         return
      end select

      maxit = 2000
      if (present(max_iterations)) maxit = max_iterations
      tol = 1.0e-10_dp
      if (present(tolerance)) tol = tolerance
      alpha = 1.0_dp
      gamma = 2.0_dp
      rho = 0.5_dp
      sigma = 0.5_dp

      allocate(simplex(n, n+1), values(n+1), centroid(n), xr(n), xe(n), xc(n))
      simplex(:, 1) = min(max(initial_model%parameters(:n), lower_bound), upper_bound)
      do i = 1, n
         simplex(:, i+1) = simplex(:, 1)
         simplex(i, i+1) = min(upper_bound(i), simplex(i, i+1) + step(i))
         if (abs(simplex(i, i+1) - simplex(i, 1)) <= tiny(1.0_dp)) then
            simplex(i, i+1) = max(lower_bound(i), simplex(i, i+1) - step(i))
         end if
      end do
      do i = 1, n+1
         values(i) = model_sse(initial_model%method, simplex(:, i), terms_years, rates)
      end do

      do iter = 1, maxit
         call order_simplex(values, best, worst, second_worst)
         simplex_spread = maxval(abs(simplex - spread(simplex(:, best), 2, n+1)))
         if (maxval(abs(values - values(best))) <= tol * (1.0_dp + abs(values(best))) .and. &
             simplex_spread <= sqrt(tol) * (1.0_dp + maxval(abs(simplex(:, best))))) exit

         centroid = (sum(simplex, dim=2) - simplex(:, worst)) / real(n, dp)
         xr = min(max(centroid + alpha * (centroid - simplex(:, worst)), lower_bound), upper_bound)
         fnew = model_sse(initial_model%method, xr, terms_years, rates)
         if (fnew < values(best)) then
            xe = min(max(centroid + gamma * (xr - centroid), lower_bound), upper_bound)
            if (model_sse(initial_model%method, xe, terms_years, rates) < fnew) then
               simplex(:, worst) = xe
               values(worst) = model_sse(initial_model%method, xe, terms_years, rates)
            else
               simplex(:, worst) = xr
               values(worst) = fnew
            end if
         else if (fnew < values(second_worst)) then
            simplex(:, worst) = xr
            values(worst) = fnew
         else
            if (fnew < values(worst)) then
               xc = min(max(centroid + rho * (xr - centroid), lower_bound), upper_bound)
            else
               xc = min(max(centroid + rho * (simplex(:, worst) - centroid), lower_bound), upper_bound)
            end if
            fnew = model_sse(initial_model%method, xc, terms_years, rates)
            if (fnew < values(worst)) then
               simplex(:, worst) = xc
               values(worst) = fnew
            else
               do i = 1, n+1
                  if (i == best) cycle
                  simplex(:, i) = min(max(simplex(:, best) + sigma * &
                     (simplex(:, i) - simplex(:, best)), lower_bound), upper_bound)
                  values(i) = model_sse(initial_model%method, simplex(:, i), terms_years, rates)
               end do
            end if
         end if
      end do

      call order_simplex(values, best, worst, second_worst)
      fit%model = initial_model
      fit%model%parameters(:n) = simplex(:, best)
      fit%objective = values(best)
      fit%iterations = min(iter, maxit)
      if (iter <= maxit) then
         fit%status = FI_OK
      else
         fit%status = FI_NO_CONVERGENCE
      end if
   end function fit_interpolation_model

   pure real(dp) function model_sse(method, parameters, terms, rates) result(value)
      integer, intent(in) :: method
      real(dp), intent(in) :: parameters(:), terms(:), rates(:)
      real(dp), allocatable :: fitted(:)
      allocate(fitted(size(rates)))
      select case (method)
      case (INTERP_NELSON_SIEGEL)
         fitted = nelson_siegel(terms, parameters(1), parameters(2), parameters(3), parameters(4))
      case (INTERP_NELSON_SIEGEL_SVENSSON)
         fitted = nelson_siegel_svensson(terms, parameters(1), parameters(2), parameters(3), &
                                          parameters(4), parameters(5), parameters(6))
      case default
         value = huge(1.0_dp)
         return
      end select
      value = sum((rates - fitted)**2)
   end function model_sse

   subroutine order_simplex(values, best, worst, second_worst)
      real(dp), intent(in) :: values(:)
      integer, intent(out) :: best, worst, second_worst
      integer :: i
      best = minloc(values, dim=1)
      worst = maxloc(values, dim=1)
      second_worst = best
      do i = 1, size(values)
         if (i == worst) cycle
         if (second_worst == worst .or. values(i) > values(second_worst)) second_worst = i
      end do
   end subroutine order_simplex

   real(dp) function linear_one(x, y, q, status) result(value)
      real(dp), intent(in) :: x(:), y(:), q
      integer, intent(out) :: status
      integer :: lo, hi, mid
      if (q < x(1) .or. q > x(size(x))) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         status = FI_OUT_OF_RANGE
         return
      end if
      if (abs(q - x(size(x))) <= 32.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(q))) then
         value = y(size(y))
         status = FI_OK
         return
      end if
      lo = 1
      hi = size(x)
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (x(mid) <= q) then
            lo = mid
         else
            hi = mid
         end if
      end do
      value = y(lo) + (q - x(lo)) * (y(lo+1) - y(lo)) / (x(lo+1) - x(lo))
      status = FI_OK
   end function linear_one

   function natural_second_derivatives(x, y, status) result(second)
      real(dp), intent(in) :: x(:), y(:)
      integer, intent(out) :: status
      real(dp), allocatable :: second(:), u(:)
      real(dp) :: sig, p
      integer :: i, n, k
      n = size(x)
      allocate(second(n), u(n))
      second = 0.0_dp
      u = 0.0_dp
      do i = 2, n-1
         sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
         p = sig * second(i-1) + 2.0_dp
         second(i) = (sig - 1.0_dp) / p
         u(i) = (6.0_dp * ((y(i+1)-y(i))/(x(i+1)-x(i)) - &
            (y(i)-y(i-1))/(x(i)-x(i-1))) / (x(i+1)-x(i-1)) - sig*u(i-1)) / p
      end do
      do k = n-1, 1, -1
         second(k) = second(k) * second(k+1) + u(k)
      end do
      status = FI_OK
   end function natural_second_derivatives

   real(dp) function natural_spline_one(x, y, second, q, status) result(value)
      real(dp), intent(in) :: x(:), y(:), second(:), q
      integer, intent(out) :: status
      integer :: lo, hi, mid
      real(dp) :: h, a, b
      if (q < x(1) .or. q > x(size(x))) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         status = FI_OUT_OF_RANGE
         return
      end if
      lo = 1
      hi = size(x)
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (x(mid) > q) then
            hi = mid
         else
            lo = mid
         end if
      end do
      h = x(hi) - x(lo)
      a = (x(hi) - q) / h
      b = (q - x(lo)) / h
      value = a*y(lo) + b*y(hi) + ((a**3-a)*second(lo) + (b**3-b)*second(hi))*h*h/6.0_dp
      status = FI_OK
   end function natural_spline_one

   function monotone_slopes(x, y) result(m)
      real(dp), intent(in) :: x(:), y(:)
      real(dp), allocatable :: m(:), delta(:), h(:)
      real(dp) :: w1, w2, a, b, tau
      integer :: i, n
      n = size(x)
      allocate(m(n), delta(n-1), h(n-1))
      h = x(2:) - x(:n-1)
      delta = (y(2:) - y(:n-1)) / h
      if (n == 2) then
         m = delta(1)
         return
      end if
      m(1) = ((2.0_dp*h(1)+h(2))*delta(1)-h(1)*delta(2))/(h(1)+h(2))
      if (m(1) * delta(1) <= 0.0_dp) m(1) = 0.0_dp
      if (delta(1) * delta(2) <= 0.0_dp .and. &
          abs(m(1)) > abs(3.0_dp*delta(1))) m(1) = 3.0_dp*delta(1)
      do i = 2, n-1
         if (delta(i-1) * delta(i) <= 0.0_dp) then
            m(i) = 0.0_dp
         else
            w1 = 2.0_dp*h(i) + h(i-1)
            w2 = h(i) + 2.0_dp*h(i-1)
            m(i) = (w1+w2)/(w1/delta(i-1)+w2/delta(i))
         end if
      end do
      m(n) = ((2.0_dp*h(n-1)+h(n-2))*delta(n-1)-h(n-1)*delta(n-2))/(h(n-1)+h(n-2))
      if (m(n) * delta(n-1) <= 0.0_dp) m(n) = 0.0_dp
      if (delta(n-1) * delta(n-2) <= 0.0_dp .and. &
          abs(m(n)) > abs(3.0_dp*delta(n-1))) m(n) = 3.0_dp*delta(n-1)
      do i = 1, n-1
         if (abs(delta(i)) <= tiny(1.0_dp)) then
            m(i) = 0.0_dp
            m(i+1) = 0.0_dp
         else
            a = m(i) / delta(i)
            b = m(i+1) / delta(i)
            if (a*a + b*b > 9.0_dp) then
               tau = 3.0_dp / sqrt(a*a + b*b)
               m(i) = tau*a*delta(i)
               m(i+1) = tau*b*delta(i)
            end if
         end if
      end do
   end function monotone_slopes

   real(dp) function hermite_one(x, y, m, q, status) result(value)
      real(dp), intent(in) :: x(:), y(:), m(:), q
      integer, intent(out) :: status
      integer :: lo, hi, mid
      real(dp) :: h, t, h00, h10, h01, h11
      if (q < x(1) .or. q > x(size(x))) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         status = FI_OUT_OF_RANGE
         return
      end if
      if (abs(q - x(size(x))) <= 32.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(q))) then
         value = y(size(y))
         status = FI_OK
         return
      end if
      lo = 1
      hi = size(x)
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (x(mid) <= q) then
            lo = mid
         else
            hi = mid
         end if
      end do
      h = x(lo+1) - x(lo)
      t = (q - x(lo)) / h
      h00 = (1.0_dp + 2.0_dp*t)*(1.0_dp-t)**2
      h10 = t*(1.0_dp-t)**2
      h01 = t*t*(3.0_dp-2.0_dp*t)
      h11 = t*t*(t-1.0_dp)
      value = h00*y(lo) + h10*h*m(lo) + h01*y(lo+1) + h11*h*m(lo+1)
      status = FI_OK
   end function hermite_one

   elemental pure real(dp) function expm1_safe(x) result(value)
      real(dp), intent(in) :: x
      if (abs(x) < 1.0e-5_dp) then
         value = x + 0.5_dp*x*x + x*x*x/6.0_dp + x**4/24.0_dp
      else
         value = exp(x) - 1.0_dp
      end if
   end function expm1_safe

end module fixedincome_interpolation
