! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from GARCHSK 0.1.0, Copyright (C) 2021 Kei Nakagawa.
module garchsk_estimation
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
   use garchsk_kinds, only : dp
   use garchsk_stats, only : sample_variance, sample_skewness, sample_kurtosis, covariance
   use garchsk_linalg, only : invert_matrix, symmetrize_matrix
   use garchsk_models, only : garchsk_negative_log_likelihood, gjrsk_negative_log_likelihood, &
      garchsk_parameters_valid, gjrsk_parameters_valid
   use garchsk_types, only : estimate_result
   implicit none
   private
   public :: garchsk_initial_parameters, gjrsk_initial_parameters
   public :: garchsk_estimate, gjrsk_estimate

contains

   function garchsk_initial_parameters(data) result(params)
      real(dp), intent(in) :: data(:)
      real(dp) :: params(10)
      real(dp) :: a1, v, sk, ku
      integer :: n
      n = size(data)
      v = max(sample_variance(data), 1.0e-8_dp)
      sk = sample_skewness(data)
      ku = max(sample_kurtosis(data), 1.0e-4_dp)
      if (n >= 3 .and. sample_variance(data(:n-1)) > tiny(1.0_dp)) then
         a1 = covariance(data(2:n), data(:n-1)) / sample_variance(data(:n-1))
      else
         a1 = 0.0_dp
      end if
      a1 = max(-0.95_dp, min(0.95_dp, a1))
      params = [a1, 0.09_dp*v, 0.01_dp, 0.90_dp, 0.29_dp*sk, 0.01_dp, 0.70_dp, &
         0.19_dp*ku, 0.01_dp, 0.80_dp]
      params(2) = max(params(2), 1.0e-8_dp)
      params(8) = max(params(8), 1.0e-8_dp)
   end function garchsk_initial_parameters

   function gjrsk_initial_parameters(data) result(params)
      real(dp), intent(in) :: data(:)
      real(dp) :: params(13)
      real(dp) :: a1, v, sk, ku
      integer :: n
      n = size(data)
      v = max(sample_variance(data), 1.0e-8_dp)
      sk = sample_skewness(data)
      ku = max(sample_kurtosis(data), 1.0e-4_dp)
      if (n >= 3 .and. sample_variance(data(:n-1)) > tiny(1.0_dp)) then
         a1 = covariance(data(2:n), data(:n-1)) / sample_variance(data(:n-1))
      else
         a1 = 0.0_dp
      end if
      a1 = max(-0.95_dp, min(0.95_dp, a1))
      params = [a1, 0.04_dp*v, 0.01_dp, 0.10_dp, 0.85_dp, 0.19_dp*sk, 0.01_dp, 0.10_dp, 0.70_dp, &
         0.14_dp*ku, 0.01_dp, 0.05_dp, 0.80_dp]
      params(2) = max(params(2), 1.0e-8_dp)
      params(10) = max(params(10), 1.0e-8_dp)
   end function gjrsk_initial_parameters

   function garchsk_estimate(data, max_iterations, tolerance, initial, include_constant) result(result)
      real(dp), intent(in) :: data(:)
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      real(dp), intent(in), optional :: initial(:)
      logical, intent(in), optional :: include_constant
      type(estimate_result) :: result
      real(dp) :: x0(10)
      integer :: maxit
      real(dp) :: tol
      logical :: add_constant

      maxit = 4000
      if (present(max_iterations)) maxit = max_iterations
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = tolerance
      add_constant = .false.
      if (present(include_constant)) add_constant = include_constant
      if (present(initial)) then
         if (size(initial) /= 10) then
            result%message = 'garchsk_estimate: initial vector must have length 10'
            return
         end if
         x0 = initial
      else
         x0 = garchsk_initial_parameters(data)
      end if
      call estimate_model(data, 1, x0, maxit, tol, add_constant, result)
   end function garchsk_estimate

   function gjrsk_estimate(data, max_iterations, tolerance, initial, include_constant) result(result)
      real(dp), intent(in) :: data(:)
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      real(dp), intent(in), optional :: initial(:)
      logical, intent(in), optional :: include_constant
      type(estimate_result) :: result
      real(dp) :: x0(13)
      integer :: maxit
      real(dp) :: tol
      logical :: add_constant

      maxit = 5000
      if (present(max_iterations)) maxit = max_iterations
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = tolerance
      add_constant = .false.
      if (present(include_constant)) add_constant = include_constant
      if (present(initial)) then
         if (size(initial) /= 13) then
            result%message = 'gjrsk_estimate: initial vector must have length 13'
            return
         end if
         x0 = initial
      else
         x0 = gjrsk_initial_parameters(data)
      end if
      call estimate_model(data, 2, x0, maxit, tol, add_constant, result)
   end function gjrsk_estimate

   subroutine estimate_model(data, model, x0, maxit, tol, add_constant, result)
      real(dp), intent(in) :: data(:), x0(:), tol
      integer, intent(in) :: model, maxit
      logical, intent(in) :: add_constant
      type(estimate_result), intent(out) :: result
      real(dp), allocatable :: best(:), hessian(:, :), covariance_matrix(:, :)
      real(dp) :: best_value, nan_value
      logical :: inverse_ok
      integer :: npar, i

      nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
      npar = size(x0)
      allocate(result%params(npar), result%standard_errors(npar), result%upstream_standard_errors(npar), &
         result%t_statistics(npar))
      result%standard_errors = nan_value
      result%upstream_standard_errors = nan_value
      result%t_statistics = nan_value
      if (size(data) < max(20, 2*npar) .or. sample_variance(data) <= tiny(1.0_dp)) then
         result%message = 'estimate_model: insufficient or constant data'
         result%params = x0
         return
      end if

      call nelder_mead(data, model, x0, maxit, tol, add_constant, best, best_value, result%iterations, &
         result%evaluations, result%converged)
      result%params = best
      result%negative_log_likelihood = best_value
      result%log_likelihood = -best_value
      result%aic = 2.0_dp * best_value + 2.0_dp * real(npar, dp)
      result%bic = 2.0_dp * best_value + real(npar, dp) * log(real(size(data), dp))
      result%upstream_aic = -2.0_dp * best_value + 2.0_dp * real(npar, dp)
      result%upstream_bic = -2.0_dp * best_value + real(npar, dp) * log(real(size(data), dp))

      call numerical_hessian(data, model, best, add_constant, hessian)
      do i = 1, npar
         if (hessian(i, i) >= 0.0_dp .and. ieee_is_finite(hessian(i, i))) then
            result%upstream_standard_errors(i) = sqrt(hessian(i, i) / real(size(data), dp))
         end if
      end do
      call invert_matrix(hessian, covariance_matrix, inverse_ok)
      if (inverse_ok) then
         call symmetrize_matrix(covariance_matrix)
         do i = 1, npar
            if (covariance_matrix(i, i) > 0.0_dp .and. ieee_is_finite(covariance_matrix(i, i))) then
               result%standard_errors(i) = sqrt(covariance_matrix(i, i))
               result%t_statistics(i) = best(i) / result%standard_errors(i)
            end if
         end do
         result%covariance_available = all(ieee_is_finite(result%standard_errors))
      end if
      if (result%converged) then
         result%message = 'ok'
      else
         result%message = 'maximum iterations reached; best feasible result returned'
      end if
   end subroutine estimate_model

   subroutine nelder_mead(data, model, x0, maxit, tol, add_constant, best, best_value, iterations, evaluations, converged)
      real(dp), intent(in) :: data(:), x0(:), tol
      integer, intent(in) :: model, maxit
      logical, intent(in) :: add_constant
      real(dp), allocatable, intent(out) :: best(:)
      real(dp), intent(out) :: best_value
      integer, intent(out) :: iterations, evaluations
      logical, intent(out) :: converged
      real(dp), allocatable :: simplex(:, :), f(:), centroid(:), xr(:), xe(:), xc(:), step(:)
      real(dp) :: fr, fe, fc, spread_x, spread_f
      integer :: n, j

      n = size(x0)
      allocate(simplex(n, n+1), f(n+1), centroid(n), xr(n), xe(n), xc(n), step(n), best(n))
      step = 0.08_dp * (abs(x0) + 0.05_dp)
      simplex(:, 1) = x0
      do j = 2, n + 1
         simplex(:, j) = x0
         simplex(j-1, j) = simplex(j-1, j) + step(j-1)
      end do
      do j = 1, n + 1
         f(j) = objective_value(simplex(:, j), data, model, add_constant)
      end do
      evaluations = n + 1
      converged = .false.

      do iterations = 1, maxit
         call sort_simplex(simplex, f)
         spread_f = maxval(abs(f(2:) - f(1)))
         spread_x = maxval(abs(simplex(:, 2:) - spread(simplex(:, 1), 2, n)))
         if (spread_f <= tol * (1.0_dp + abs(f(1))) .and. &
            spread_x <= sqrt(tol) * (1.0_dp + maxval(abs(simplex(:, 1))))) then
            converged = .true.
            exit
         end if
         centroid = sum(simplex(:, :n), dim=2) / real(n, dp)
         xr = centroid + (centroid - simplex(:, n+1))
         fr = objective_value(xr, data, model, add_constant)
         evaluations = evaluations + 1
         if (fr < f(1)) then
            xe = centroid + 2.0_dp * (xr - centroid)
            fe = objective_value(xe, data, model, add_constant)
            evaluations = evaluations + 1
            if (fe < fr) then
               simplex(:, n+1) = xe
               f(n+1) = fe
            else
               simplex(:, n+1) = xr
               f(n+1) = fr
            end if
         else if (fr < f(n)) then
            simplex(:, n+1) = xr
            f(n+1) = fr
         else
            if (fr < f(n+1)) then
               xc = centroid + 0.5_dp * (xr - centroid)
            else
               xc = centroid + 0.5_dp * (simplex(:, n+1) - centroid)
            end if
            fc = objective_value(xc, data, model, add_constant)
            evaluations = evaluations + 1
            if (fc < min(fr, f(n+1))) then
               simplex(:, n+1) = xc
               f(n+1) = fc
            else
               do j = 2, n + 1
                  simplex(:, j) = simplex(:, 1) + 0.5_dp * (simplex(:, j) - simplex(:, 1))
                  f(j) = objective_value(simplex(:, j), data, model, add_constant)
               end do
               evaluations = evaluations + n
            end if
         end if
      end do
      call sort_simplex(simplex, f)
      best = simplex(:, 1)
      best_value = f(1)
   end subroutine nelder_mead

   real(dp) function objective_value(x, data, model, add_constant) result(value)
      real(dp), intent(in) :: x(:), data(:)
      integer, intent(in) :: model
      logical, intent(in) :: add_constant
      if (model == 1) then
         value = garchsk_negative_log_likelihood(x, data, add_constant)
      else
         value = gjrsk_negative_log_likelihood(x, data, add_constant)
      end if
   end function objective_value

   subroutine sort_simplex(simplex, f)
      real(dp), intent(inout) :: simplex(:, :), f(:)
      real(dp), allocatable :: temp(:)
      real(dp) :: tf
      integer :: i, j, k
      allocate(temp(size(simplex, 1)))
      do i = 1, size(f) - 1
         k = i
         do j = i + 1, size(f)
            if (f(j) < f(k)) k = j
         end do
         if (k /= i) then
            tf = f(i)
            f(i) = f(k)
            f(k) = tf
            temp = simplex(:, i)
            simplex(:, i) = simplex(:, k)
            simplex(:, k) = temp
         end if
      end do
   end subroutine sort_simplex

   subroutine numerical_hessian(data, model, x, add_constant, hessian)
      real(dp), intent(in) :: data(:), x(:)
      integer, intent(in) :: model
      logical, intent(in) :: add_constant
      real(dp), allocatable, intent(out) :: hessian(:, :)
      real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:), step(:)
      real(dp) :: f0
      integer :: n, i, j

      n = size(x)
      allocate(hessian(n, n), xp(n), xm(n), xpp(n), xpm(n), xmp(n), xmm(n), step(n))
      step = epsilon(1.0_dp)**0.25_dp * (abs(x) + 1.0_dp)
      f0 = objective_value(x, data, model, add_constant)
      do i = 1, n
         xp = x
         xm = x
         xp(i) = xp(i) + step(i)
         xm(i) = xm(i) - step(i)
         hessian(i, i) = (objective_value(xp, data, model, add_constant) - 2.0_dp*f0 + &
            objective_value(xm, data, model, add_constant)) / step(i)**2
         do j = i + 1, n
            xpp = x
            xpm = x
            xmp = x
            xmm = x
            xpp(i) = xpp(i) + step(i)
            xpp(j) = xpp(j) + step(j)
            xpm(i) = xpm(i) + step(i)
            xpm(j) = xpm(j) - step(j)
            xmp(i) = xmp(i) - step(i)
            xmp(j) = xmp(j) + step(j)
            xmm(i) = xmm(i) - step(i)
            xmm(j) = xmm(j) - step(j)
            hessian(i, j) = (objective_value(xpp, data, model, add_constant) - &
               objective_value(xpm, data, model, add_constant) - &
               objective_value(xmp, data, model, add_constant) + &
               objective_value(xmm, data, model, add_constant)) / (4.0_dp*step(i)*step(j))
            hessian(j, i) = hessian(i, j)
         end do
      end do
   end subroutine numerical_hessian

end module garchsk_estimation
