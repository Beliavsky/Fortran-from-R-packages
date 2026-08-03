! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module sandwich_hac
   use sandwich_kinds, only : dp
   use sandwich_status, only : SANDWICH_SUCCESS, SANDWICH_INVALID_ARGUMENT, &
      SANDWICH_DIMENSION_MISMATCH, SANDWICH_NUMERICAL_FAILURE
   use sandwich_utils, only : lowercase
   use sandwich_linalg, only : inverse_matrix, identity_matrix
   use sandwich_kernels, only : kernel_weights
   use sandwich_auxiliary, only : isoacf
   use sandwich_core, only : sandwich_covariance
   implicit none
   private

   type, public :: hac_diagnostics
      real(dp) :: bias_correction = 1.0_dp
      real(dp) :: degrees_freedom = 0.0_dp
      real(dp) :: bandwidth = 0.0_dp
      integer :: effective_n = 0
      integer :: nweights = 0
      integer :: prewhite_order = 0
   end type hac_diagnostics

   public :: meat_hac, vcov_hac
   public :: bandwidth_andrews, andrews_weights
   public :: bandwidth_newey_west, newey_west_weights
   public :: lumley_weights, long_run_variance
   public :: prewhite_var

contains

   subroutine prewhite_var(scores, order, residuals, recolor, status)
      real(dp), intent(in) :: scores(:, :)
      integer, intent(in) :: order
      real(dp), allocatable, intent(out) :: residuals(:, :), recolor(:, :)
      integer, intent(out), optional :: status
      real(dp), allocatable :: z(:, :), y(:, :), ztz(:, :), zty(:, :), coef(:, :)
      real(dp), allocatable :: sum_a(:, :), eye(:, :)
      integer :: n, k, m, t, lag, info

      n = size(scores, 1)
      k = size(scores, 2)
      if (order < 0 .or. n <= order .or. k <= 0) then
         allocate(residuals(0, 0), recolor(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      if (order == 0) then
         allocate(residuals(n, k))
         residuals = scores
         recolor = identity_matrix(k)
         if (present(status)) status = SANDWICH_SUCCESS
         return
      end if

      m = n - order
      allocate(z(m, k * order), y(m, k))
      do t = 1, m
         y(t, :) = scores(order + t, :)
         do lag = 1, order
            z(t, (lag - 1) * k + 1:lag * k) = scores(order + t - lag, :)
         end do
      end do
      ztz = matmul(transpose(z), z)
      zty = matmul(transpose(z), y)
      call solve_multiple(ztz, zty, coef, info)
      if (info /= SANDWICH_SUCCESS) then
         allocate(residuals(0, 0), recolor(0, 0))
         if (present(status)) status = info
         return
      end if

      residuals = y - matmul(z, coef)
      allocate(sum_a(k, k))
      sum_a = 0.0_dp
      do lag = 1, order
         sum_a = sum_a + transpose(coef((lag - 1) * k + 1:lag * k, :))
      end do
      eye = identity_matrix(k)
      call inverse_matrix(eye - sum_a, recolor, info)
      if (info /= SANDWICH_SUCCESS) then
         deallocate(residuals)
         allocate(residuals(0, 0), recolor(0, 0))
      end if
      if (present(status)) status = info
   contains
      subroutine solve_multiple(a, b, x, solve_status)
         use sandwich_linalg, only : solve_linear
         real(dp), intent(in) :: a(:, :), b(:, :)
         real(dp), allocatable, intent(out) :: x(:, :)
         integer, intent(out) :: solve_status
         call solve_linear(a, b, x, solve_status)
      end subroutine solve_multiple
   end subroutine prewhite_var

   subroutine meat_hac(scores, weights, meat_matrix, status, adjust, prewhite_order, diagnostics)
      real(dp), intent(in) :: scores(:, :), weights(:)
      real(dp), allocatable, intent(out) :: meat_matrix(:, :)
      integer, intent(out), optional :: status
      logical, intent(in), optional :: adjust
      integer, intent(in), optional :: prewhite_order
      type(hac_diagnostics), intent(out), optional :: diagnostics
      real(dp), allocatable :: u(:, :), d(:, :), utu(:, :)
      real(dp) :: wsum, w2sum
      integer :: n_orig, n, k, p, nw, lag, info
      logical :: use_adjust

      n_orig = size(scores, 1)
      k = size(scores, 2)
      p = 0
      if (present(prewhite_order)) p = prewhite_order
      use_adjust = .true.
      if (present(adjust)) use_adjust = adjust

      if (n_orig <= 0 .or. k <= 0 .or. p < 0 .or. n_orig <= p .or. size(weights) <= 0) then
         allocate(meat_matrix(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      if (use_adjust .and. n_orig <= k) then
         allocate(meat_matrix(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if

      call prewhite_var(scores, p, u, d, info)
      if (info /= SANDWICH_SUCCESS) then
         allocate(meat_matrix(0, 0))
         if (present(status)) status = info
         return
      end if
      n = size(u, 1)
      nw = min(size(weights), n)

      utu = 0.5_dp * weights(1) * matmul(transpose(u), u)
      wsum = real(n, dp) * weights(1) / 2.0_dp
      w2sum = real(n, dp) * weights(1)**2 / 2.0_dp
      do lag = 1, nw - 1
         utu = utu + weights(lag + 1) * matmul(transpose(u(1:n - lag, :)), u(1 + lag:n, :))
         wsum = wsum + real(n - lag, dp) * weights(lag + 1)
         w2sum = w2sum + real(n - lag, dp) * weights(lag + 1)**2
      end do
      utu = utu + transpose(utu)

      if (use_adjust) utu = real(n_orig, dp) / real(n_orig - k, dp) * utu
      if (p > 0) utu = matmul(d, matmul(utu, transpose(d)))
      meat_matrix = utu / real(n_orig, dp)
      meat_matrix = 0.5_dp * (meat_matrix + transpose(meat_matrix))

      if (present(diagnostics)) then
         wsum = 2.0_dp * wsum
         w2sum = 2.0_dp * w2sum
         diagnostics%bias_correction = real(n, dp)**2 / (real(n, dp)**2 - wsum)
         diagnostics%degrees_freedom = real(n, dp)**2 / w2sum
         diagnostics%effective_n = n
         diagnostics%nweights = nw
         diagnostics%prewhite_order = p
      end if
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine meat_hac

   subroutine vcov_hac(scores, bread, weights, covariance, status, adjust, prewhite_order, diagnostics)
      real(dp), intent(in) :: scores(:, :), bread(:, :), weights(:)
      real(dp), allocatable, intent(out) :: covariance(:, :)
      integer, intent(out), optional :: status
      logical, intent(in), optional :: adjust
      integer, intent(in), optional :: prewhite_order
      type(hac_diagnostics), intent(out), optional :: diagnostics
      real(dp), allocatable :: meat_matrix(:, :)
      type(hac_diagnostics) :: diag_local
      integer :: info, p
      logical :: use_adjust

      p = 0
      if (present(prewhite_order)) p = prewhite_order
      use_adjust = .true.
      if (present(adjust)) use_adjust = adjust
      call meat_hac(scores, weights, meat_matrix, info, use_adjust, p, diag_local)
      if (info /= SANDWICH_SUCCESS) then
         allocate(covariance(0, 0))
         if (present(status)) status = info
         return
      end if
      call sandwich_covariance(bread, meat_matrix, size(scores, 1), covariance, info)
      if (present(diagnostics)) diagnostics = diag_local
      if (present(status)) status = info
   end subroutine vcov_hac

   subroutine fit_ar1(x, rho, sigma, status)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: rho, sigma
      integer, intent(out) :: status
      real(dp) :: denominator
      integer :: n

      n = size(x)
      if (n < 3) then
         rho = 0.0_dp
         sigma = 0.0_dp
         status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      denominator = dot_product(x(1:n - 1), x(1:n - 1))
      if (denominator <= tiny(1.0_dp)) then
         rho = 0.0_dp
         sigma = 0.0_dp
         status = SANDWICH_NUMERICAL_FAILURE
         return
      end if
      rho = dot_product(x(2:n), x(1:n - 1)) / denominator
      rho = max(-0.999_dp, min(0.999_dp, rho))
      sigma = sqrt(sum((x(2:n) - rho * x(1:n - 1))**2) / real(n - 1, dp))
      status = SANDWICH_SUCCESS
   end subroutine fit_ar1

   real(dp) function arma11_objective(x, rho, psi) result(value)
      real(dp), intent(in) :: x(:), rho, psi
      real(dp) :: previous_error, current_error
      integer :: t

      if (abs(rho) >= 0.999_dp .or. abs(psi) >= 0.999_dp) then
         value = huge(1.0_dp) / 1000.0_dp
         return
      end if
      previous_error = x(1)
      value = 0.0_dp
      do t = 2, size(x)
         current_error = x(t) - rho * x(t - 1) - psi * previous_error
         value = value + current_error * current_error
         previous_error = current_error
      end do
   end function arma11_objective

   subroutine fit_arma11_css(x, rho, psi, sigma, status)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: rho, psi, sigma
      integer, intent(out) :: status
      real(dp) :: simplex(2, 3), f(3), centroid(2), reflected(2), expanded(2), contracted(2)
      real(dp) :: fr, fe, fc, spread
      integer :: i, j, iter, best, worst, second

      if (size(x) < 6) then
         rho = 0.0_dp
         psi = 0.0_dp
         sigma = 0.0_dp
         status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      call fit_ar1(x, rho, sigma, status)
      if (status /= SANDWICH_SUCCESS) return
      psi = 0.0_dp
      simplex(:, 1) = [rho, psi]
      simplex(:, 2) = [max(-0.95_dp, min(0.95_dp, rho + 0.08_dp)), psi]
      simplex(:, 3) = [rho, 0.08_dp]
      do i = 1, 3
         f(i) = arma11_objective(x, simplex(1, i), simplex(2, i))
      end do

      do iter = 1, 300
         call order_three(f, best, second, worst)
         spread = maxval(abs(simplex(:, worst) - simplex(:, best)))
         if (spread < 1.0e-8_dp .and. abs(f(worst) - f(best)) < 1.0e-10_dp * max(1.0_dp, f(best))) exit
         centroid = 0.5_dp * (simplex(:, best) + simplex(:, second))
         reflected = centroid + (centroid - simplex(:, worst))
         reflected = max(-0.998_dp, min(0.998_dp, reflected))
         fr = arma11_objective(x, reflected(1), reflected(2))
         if (fr < f(best)) then
            expanded = centroid + 2.0_dp * (reflected - centroid)
            expanded = max(-0.998_dp, min(0.998_dp, expanded))
            fe = arma11_objective(x, expanded(1), expanded(2))
            if (fe < fr) then
               simplex(:, worst) = expanded
               f(worst) = fe
            else
               simplex(:, worst) = reflected
               f(worst) = fr
            end if
         else if (fr < f(second)) then
            simplex(:, worst) = reflected
            f(worst) = fr
         else
            if (fr < f(worst)) then
               contracted = centroid + 0.5_dp * (reflected - centroid)
            else
               contracted = centroid + 0.5_dp * (simplex(:, worst) - centroid)
            end if
            contracted = max(-0.998_dp, min(0.998_dp, contracted))
            fc = arma11_objective(x, contracted(1), contracted(2))
            if (fc < min(fr, f(worst))) then
               simplex(:, worst) = contracted
               f(worst) = fc
            else
               do j = 1, 3
                  if (j /= best) then
                     simplex(:, j) = simplex(:, best) + 0.5_dp * (simplex(:, j) - simplex(:, best))
                     f(j) = arma11_objective(x, simplex(1, j), simplex(2, j))
                  end if
               end do
            end if
         end if
      end do
      call order_three(f, best, second, worst)
      rho = simplex(1, best)
      psi = simplex(2, best)
      sigma = sqrt(f(best) / real(size(x) - 1, dp))
      status = SANDWICH_SUCCESS
   contains
      subroutine order_three(values, first, middle, last)
         real(dp), intent(in) :: values(3)
         integer, intent(out) :: first, middle, last
         integer :: idx(3), a, b, tmp_index
         idx = [1, 2, 3]
         do a = 1, 2
            do b = a + 1, 3
               if (values(idx(b)) < values(idx(a))) then
                  tmp_index = idx(a)
                  idx(a) = idx(b)
                  idx(b) = tmp_index
               end if
            end do
         end do
         first = idx(1)
         middle = idx(2)
         last = idx(3)
      end subroutine order_three
   end subroutine fit_arma11_css

   subroutine bandwidth_andrews(scores, kernel, bandwidth, status, approx, series_weights, prewhite_order)
      real(dp), intent(in) :: scores(:, :)
      character(len=*), intent(in) :: kernel
      real(dp), intent(out) :: bandwidth
      integer, intent(out), optional :: status
      character(len=*), intent(in), optional :: approx
      real(dp), intent(in), optional :: series_weights(:)
      integer, intent(in), optional :: prewhite_order
      real(dp), allocatable :: u(:, :), d(:, :), w(:)
      real(dp) :: denum, alpha1, alpha2, rho, psi, sigma, numerator1, numerator2
      integer :: n, k, j, p, info
      character(len=:), allocatable :: approximation, kind

      n = size(scores, 1)
      k = size(scores, 2)
      p = 1
      if (present(prewhite_order)) p = prewhite_order
      approximation = 'ar(1)'
      if (present(approx)) approximation = trim(lowercase(approx))
      kind = trim(lowercase(kernel))

      if (n <= p + 2 .or. k <= 0) then
         bandwidth = 0.0_dp
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      allocate(w(k))
      w = 1.0_dp
      if (present(series_weights)) then
         if (size(series_weights) /= k) then
            bandwidth = 0.0_dp
            if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
            return
         end if
         w = series_weights
      end if
      if (maxval(abs(w)) <= tiny(1.0_dp)) w = 1.0_dp

      call prewhite_var(scores, p, u, d, info)
      if (info /= SANDWICH_SUCCESS) then
         bandwidth = 0.0_dp
         if (present(status)) status = info
         return
      end if
      n = size(u, 1)
      denum = 0.0_dp
      numerator1 = 0.0_dp
      numerator2 = 0.0_dp
      do j = 1, k
         if (approximation == 'ar(1)' .or. approximation == 'ar1') then
            call fit_ar1(u(:, j), rho, sigma, info)
            if (info /= SANDWICH_SUCCESS) then
               bandwidth = 0.0_dp
               if (present(status)) status = info
               return
            end if
            denum = denum + w(j) * (sigma / (1.0_dp - rho))**4
            numerator2 = numerator2 + w(j) * 4.0_dp * rho**2 * sigma**4 / (1.0_dp - rho)**8
            numerator1 = numerator1 + w(j) * 4.0_dp * rho**2 * sigma**4 / &
               ((1.0_dp - rho)**6 * (1.0_dp + rho)**2)
         else if (approximation == 'arma(1,1)' .or. approximation == 'arma11') then
            call fit_arma11_css(u(:, j), rho, psi, sigma, info)
            if (info /= SANDWICH_SUCCESS) then
               bandwidth = 0.0_dp
               if (present(status)) status = info
               return
            end if
            denum = denum + w(j) * ((1.0_dp + psi) * sigma / (1.0_dp - rho))**4
            numerator2 = numerator2 + w(j) * 4.0_dp * ((1.0_dp + rho * psi) * &
               (rho + psi))**2 * sigma**4 / (1.0_dp - rho)**8
            numerator1 = numerator1 + w(j) * 4.0_dp * ((1.0_dp + rho * psi) * &
               (rho + psi))**2 * sigma**4 / ((1.0_dp - rho)**6 * (1.0_dp + rho)**2)
         else
            bandwidth = 0.0_dp
            if (present(status)) status = SANDWICH_INVALID_ARGUMENT
            return
         end if
      end do
      if (denum <= tiny(1.0_dp)) then
         bandwidth = 0.0_dp
         if (present(status)) status = SANDWICH_NUMERICAL_FAILURE
         return
      end if
      alpha1 = numerator1 / denum
      alpha2 = numerator2 / denum

      select case (kind)
      case ('truncated', 'uniform')
         bandwidth = 0.6611_dp * (real(n, dp) * alpha2)**0.2_dp
      case ('bartlett', 'triangular')
         bandwidth = 1.1447_dp * (real(n, dp) * alpha1)**(1.0_dp / 3.0_dp)
      case ('parzen')
         bandwidth = 2.6614_dp * (real(n, dp) * alpha2)**0.2_dp
      case ('tukey-hanning', 'tukey hanning', 'tukey_hanning')
         bandwidth = 1.7462_dp * (real(n, dp) * alpha2)**0.2_dp
      case ('quadratic spectral', 'quadratic-spectral', 'quadratic_spectral', 'qs')
         bandwidth = 1.3221_dp * (real(n, dp) * alpha2)**0.2_dp
      case default
         bandwidth = 0.0_dp
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end select
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine bandwidth_andrews

   subroutine andrews_weights(scores, kernel, weights, status, bandwidth, approx, &
      series_weights, prewhite_order, tolerance)
      real(dp), intent(in) :: scores(:, :)
      character(len=*), intent(in) :: kernel
      real(dp), allocatable, intent(out) :: weights(:)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: bandwidth
      character(len=*), intent(in), optional :: approx
      real(dp), intent(in), optional :: series_weights(:)
      integer, intent(in), optional :: prewhite_order
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: x(:), all_weights(:)
      real(dp) :: bw, tol
      integer :: p, n, i, last, info
      character(len=32) :: approximation

      p = 1
      if (present(prewhite_order)) p = prewhite_order
      approximation = 'AR(1)'
      if (present(approx)) approximation = approx
      tol = 1.0e-7_dp
      if (present(tolerance)) tol = tolerance
      if (present(bandwidth)) then
         bw = bandwidth
      else if (present(series_weights)) then
         call bandwidth_andrews(scores, kernel, bw, info, approximation, series_weights, p)
      else
         call bandwidth_andrews(scores, kernel, bw, info, approximation, prewhite_order = p)
      end if
      if (.not. present(bandwidth)) then
         if (info /= SANDWICH_SUCCESS) then
            allocate(weights(0))
            if (present(status)) status = info
            return
         end if
      end if
      if (bw <= 0.0_dp) then
         allocate(weights(0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if

      n = size(scores, 1) - p
      allocate(x(n))
      do i = 1, n
         x(i) = real(i - 1, dp) / bw
      end do
      call kernel_weights(x, kernel, all_weights, info)
      if (info /= SANDWICH_SUCCESS) then
         allocate(weights(0))
         if (present(status)) status = info
         return
      end if
      last = 1
      do i = 1, n
         if (abs(all_weights(i)) > tol) last = i
      end do
      allocate(weights(last))
      weights = all_weights(1:last)
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine andrews_weights

   subroutine bandwidth_newey_west(scores, kernel, bandwidth, status, series_weights, prewhite_order)
      real(dp), intent(in) :: scores(:, :)
      character(len=*), intent(in) :: kernel
      real(dp), intent(out) :: bandwidth
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: series_weights(:)
      integer, intent(in), optional :: prewhite_order
      real(dp), allocatable :: u(:, :), d(:, :), w(:), hw(:)
      real(dp) :: mrate, qrate, s0, s1, s2, ratio, sigma_lag
      integer :: n_orig, n, k, p, m, lag, info
      character(len=:), allocatable :: kind

      n_orig = size(scores, 1)
      k = size(scores, 2)
      p = 1
      if (present(prewhite_order)) p = prewhite_order
      kind = trim(lowercase(kernel))
      select case (kind)
      case ('bartlett', 'triangular')
         mrate = 2.0_dp / 9.0_dp
      case ('parzen')
         mrate = 4.0_dp / 25.0_dp
      case ('quadratic spectral', 'quadratic-spectral', 'quadratic_spectral', 'qs')
         mrate = 2.0_dp / 25.0_dp
      case default
         bandwidth = 0.0_dp
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end select
      if (n_orig <= p + 2 .or. k <= 0) then
         bandwidth = 0.0_dp
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if

      allocate(w(k))
      w = 1.0_dp
      if (present(series_weights)) then
         if (size(series_weights) /= k) then
            bandwidth = 0.0_dp
            if (present(status)) status = SANDWICH_DIMENSION_MISMATCH
            return
         end if
         w = series_weights
      end if
      if (all(w <= 0.0_dp)) w = 1.0_dp

      m = floor(merge(3.0_dp, 4.0_dp, p > 0) * (real(n_orig, dp) / 100.0_dp)**mrate)
      call prewhite_var(scores, p, u, d, info)
      if (info /= SANDWICH_SUCCESS) then
         bandwidth = 0.0_dp
         if (present(status)) status = info
         return
      end if
      n = size(u, 1)
      m = max(0, min(m, n - 1))
      allocate(hw(n))
      hw = matmul(u, w)
      s0 = dot_product(hw, hw) / real(n, dp)
      s1 = 0.0_dp
      s2 = 0.0_dp
      do lag = 1, m
         sigma_lag = dot_product(hw(1:n - lag), hw(1 + lag:n)) / real(n, dp)
         s0 = s0 + 2.0_dp * sigma_lag
         s1 = s1 + 2.0_dp * real(lag, dp) * sigma_lag
         s2 = s2 + 2.0_dp * real(lag * lag, dp) * sigma_lag
      end do
      if (abs(s0) <= tiny(1.0_dp)) then
         bandwidth = 0.0_dp
         if (present(status)) status = SANDWICH_NUMERICAL_FAILURE
         return
      end if

      if (kind == 'bartlett' .or. kind == 'triangular') then
         qrate = 1.0_dp / 3.0_dp
         ratio = (s1 / s0)**2
         bandwidth = 1.1447_dp * max(ratio, 0.0_dp)**qrate
      else if (kind == 'parzen') then
         qrate = 0.2_dp
         ratio = (s2 / s0)**2
         bandwidth = 2.6614_dp * max(ratio, 0.0_dp)**qrate
      else
         qrate = 0.2_dp
         ratio = (s2 / s0)**2
         bandwidth = 1.3221_dp * max(ratio, 0.0_dp)**qrate
      end if
      bandwidth = bandwidth * real(n + p, dp)**qrate
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine bandwidth_newey_west

   subroutine newey_west_weights(scores, weights, status, lag, prewhite_order, series_weights)
      real(dp), intent(in) :: scores(:, :)
      real(dp), allocatable, intent(out) :: weights(:)
      integer, intent(out), optional :: status
      integer, intent(in), optional :: lag, prewhite_order
      real(dp), intent(in), optional :: series_weights(:)
      real(dp) :: bw
      integer :: selected_lag, p, j, info

      p = 1
      if (present(prewhite_order)) p = prewhite_order
      if (present(lag)) then
         selected_lag = lag
      else if (present(series_weights)) then
         call bandwidth_newey_west(scores, 'Bartlett', bw, info, series_weights, p)
         selected_lag = floor(bw)
      else
         call bandwidth_newey_west(scores, 'Bartlett', bw, info, prewhite_order = p)
         selected_lag = floor(bw)
      end if
      if (.not. present(lag)) then
         if (info /= SANDWICH_SUCCESS) then
            allocate(weights(0))
            if (present(status)) status = info
            return
         end if
      end if
      if (selected_lag < 0 .or. selected_lag >= size(scores, 1) - p) then
         allocate(weights(0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      allocate(weights(selected_lag + 1))
      do j = 0, selected_lag
         weights(j + 1) = 1.0_dp - real(j, dp) / real(selected_lag + 1, dp)
      end do
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine newey_west_weights

   subroutine lumley_weights(residuals, weights, status, method, c_value, tolerance, weave1)
      real(dp), intent(in) :: residuals(:)
      real(dp), allocatable, intent(out) :: weights(:)
      integer, intent(out), optional :: status
      character(len=*), intent(in), optional :: method
      real(dp), intent(in), optional :: c_value, tolerance
      logical, intent(in), optional :: weave1
      real(dp), allocatable :: rho(:), all_weights(:)
      real(dp) :: c, tol
      integer :: n, lag, last, info
      character(len=16) :: selected_method
      logical :: use_weave1

      n = size(residuals)
      selected_method = 'truncate'
      if (present(method)) selected_method = trim(lowercase(method))
      tol = 1.0e-7_dp
      if (present(tolerance)) tol = tolerance
      use_weave1 = .false.
      if (present(weave1)) use_weave1 = weave1
      call isoacf(residuals, rho, info, weave1 = use_weave1)
      if (info /= SANDWICH_SUCCESS) then
         allocate(weights(0))
         if (present(status)) status = info
         return
      end if

      if (trim(selected_method) == 'truncate') then
         c = 4.0_dp
         if (present(c_value)) c = c_value
         last = 1
         do lag = 1, size(rho)
            if (rho(lag)**2 * real(n, dp) > c) last = lag
         end do
         allocate(weights(last))
         weights = 1.0_dp
      else if (trim(selected_method) == 'smooth') then
         c = 1.0_dp
         if (present(c_value)) c = c_value
         allocate(all_weights(size(rho)))
         all_weights = min(1.0_dp, c * real(n, dp) * rho**2)
         last = 1
         do lag = 1, size(all_weights)
            if (abs(all_weights(lag)) > tol) last = lag
         end do
         allocate(weights(last))
         weights = all_weights(1:last)
      else
         allocate(weights(0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine lumley_weights

   subroutine long_run_variance(series, covariance, status, estimator, prewhite_order, adjust, lag, kernel)
      real(dp), intent(in) :: series(:, :)
      real(dp), allocatable, intent(out) :: covariance(:, :)
      integer, intent(out), optional :: status
      character(len=*), intent(in), optional :: estimator
      integer, intent(in), optional :: prewhite_order, lag
      logical, intent(in), optional :: adjust
      character(len=*), intent(in), optional :: kernel
      real(dp), allocatable :: scores(:, :), weights(:), meat_matrix(:, :)
      integer :: n, k, i, p, info
      logical :: use_adjust
      character(len=32) :: selected_estimator, selected_kernel

      n = size(series, 1)
      k = size(series, 2)
      if (n <= 2 .or. k <= 0) then
         allocate(covariance(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      p = 1
      if (present(prewhite_order)) p = prewhite_order
      use_adjust = .true.
      if (present(adjust)) use_adjust = adjust
      selected_estimator = 'andrews'
      if (present(estimator)) selected_estimator = trim(lowercase(estimator))
      selected_kernel = 'Quadratic Spectral'
      if (present(kernel)) selected_kernel = kernel

      allocate(scores(n, k))
      do i = 1, k
         scores(:, i) = series(:, i) - sum(series(:, i)) / real(n, dp)
      end do
      if (trim(selected_estimator) == 'newey-west' .or. trim(selected_estimator) == 'neweywest') then
         if (present(lag)) then
            call newey_west_weights(scores, weights, info, lag = lag, prewhite_order = p)
         else
            call newey_west_weights(scores, weights, info, prewhite_order = p)
         end if
      else if (trim(selected_estimator) == 'andrews') then
         call andrews_weights(scores, selected_kernel, weights, info, prewhite_order = p)
      else
         allocate(covariance(0, 0))
         if (present(status)) status = SANDWICH_INVALID_ARGUMENT
         return
      end if
      if (info /= SANDWICH_SUCCESS) then
         allocate(covariance(0, 0))
         if (present(status)) status = info
         return
      end if
      call meat_hac(scores, weights, meat_matrix, info, use_adjust, p)
      if (info /= SANDWICH_SUCCESS) then
         allocate(covariance(0, 0))
         if (present(status)) status = info
         return
      end if
      covariance = meat_matrix / real(n, dp)
      covariance = 0.5_dp * (covariance + transpose(covariance))
      if (present(status)) status = SANDWICH_SUCCESS
   end subroutine long_run_variance

end module sandwich_hac
