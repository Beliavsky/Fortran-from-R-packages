! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! Based on fBonds, copyright its original authors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License, version 2 or later.
module fbonds_term_structure
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use fbonds_kinds, only : dp
   use fbonds_linalg, only : least_squares
   implicit none
   private
   public :: term_structure_fit
   public :: nelson_siegel_curve, svensson_curve
   public :: fit_nelson_siegel, fit_svensson

   integer, parameter :: model_ns = 1, model_svensson = 2

   type :: term_structure_fit
      character(len=24) :: model = ""
      character(len=8) :: objective = ""
      real(dp), allocatable :: parameters(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: objective_value = huge(1.0_dp)
      real(dp) :: sse = huge(1.0_dp)
      real(dp) :: mae = huge(1.0_dp)
      real(dp) :: rmse = huge(1.0_dp)
      logical :: converged = .false.
      integer :: iterations = 0
      integer :: evaluations = 0
      integer :: status = 0
   end type term_structure_fit
contains
   pure elemental function ns_loading1(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      if (abs(x) < 1.0e-6_dp) then
         value = 1.0_dp - 0.5_dp * x + x * x / 6.0_dp - x**3 / 24.0_dp
      else
         value = (1.0_dp - exp(-x)) / x
      end if
   end function ns_loading1

   pure elemental function ns_loading2(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      value = ns_loading1(x) - exp(-x)
   end function ns_loading2

   pure function nelson_siegel_curve(maturity, parameters) result(rate)
      real(dp), intent(in) :: maturity(:), parameters(4)
      real(dp) :: rate(size(maturity)), x(size(maturity))
      if (parameters(4) <= 0.0_dp) then
         rate = huge(1.0_dp)
         return
      end if
      x = maturity / parameters(4)
      rate = parameters(1) + parameters(2) * ns_loading1(x) + &
         parameters(3) * ns_loading2(x)
   end function nelson_siegel_curve

   pure function svensson_curve(maturity, parameters) result(rate)
      real(dp), intent(in) :: maturity(:), parameters(6)
      real(dp) :: rate(size(maturity)), x1(size(maturity)), x2(size(maturity))
      if (parameters(5) <= 0.0_dp .or. parameters(6) <= 0.0_dp) then
         rate = huge(1.0_dp)
         return
      end if
      x1 = maturity / parameters(5)
      x2 = maturity / parameters(6)
      rate = parameters(1) + parameters(2) * exp(-x1) + &
         parameters(3) * x1 * exp(-x1) + parameters(4) * x2 * exp(-x2)
   end function svensson_curve

   subroutine fit_nelson_siegel(rate, maturity, fit, tau_lower, tau_upper, &
      max_iterations, tolerance)
      real(dp), intent(in) :: rate(:), maturity(:)
      type(term_structure_fit), intent(out) :: fit
      real(dp), intent(in), optional :: tau_lower, tau_upper, tolerance
      integer, intent(in), optional :: max_iterations
      real(dp), allocatable :: start(:), step(:), best(:)
      real(dp) :: lower_tau, upper_tau, fbest, tol
      integer :: maxit, evals

      call initialize_fit(fit, "Nelson-Siegel", "sse", size(rate), 4)
      if (.not. valid_data(rate, maturity, 4)) then
         fit%status = -1
         return
      end if
      call choose_tau_bounds(maturity, tau_lower, tau_upper, lower_tau, upper_tau)
      allocate(start(4), step(4), best(4))
      call ns_grid_start(rate, maturity, lower_tau, upper_tau, start)
      start(4) = log(start(4))
      step(1:3) = max(0.05_dp * max(1.0_dp, abs(start(1:3))), 1.0e-4_dp)
      step(4) = 0.15_dp
      maxit = 2000
      if (present(max_iterations)) maxit = max(1, max_iterations)
      tol = 1.0e-10_dp
      if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
      call optimize_term_model(model_ns, "sse", rate, maturity, start, step, &
         log(lower_tau), log(upper_tau), best, fbest, fit%converged, &
         fit%iterations, evals, maxit, tol)
      fit%evaluations = evals
      fit%parameters(1:3) = best(1:3)
      fit%parameters(4) = exp(best(4))
      fit%objective_value = fbest
      fit%fitted = nelson_siegel_curve(maturity, fit%parameters)
      call finish_fit(rate, fit)
   end subroutine fit_nelson_siegel

   subroutine fit_svensson(rate, maturity, fit, tau_lower, tau_upper, objective, &
      max_iterations, tolerance)
      real(dp), intent(in) :: rate(:), maturity(:)
      type(term_structure_fit), intent(out) :: fit
      real(dp), intent(in), optional :: tau_lower, tau_upper, tolerance
      character(len=*), intent(in), optional :: objective
      integer, intent(in), optional :: max_iterations
      real(dp), allocatable :: start(:), step(:), best(:)
      real(dp) :: lower_tau, upper_tau, fbest, tol
      character(len=8) :: objective_name
      integer :: maxit, evals

      objective_name = "l1"
      if (present(objective)) objective_name = lower_string(trim(objective))
      if (objective_name /= "l1" .and. objective_name /= "sse") objective_name = "l1"
      call initialize_fit(fit, "Svensson", objective_name, size(rate), 6)
      if (.not. valid_data(rate, maturity, 6)) then
         fit%status = -1
         return
      end if
      call choose_tau_bounds(maturity, tau_lower, tau_upper, lower_tau, upper_tau)
      allocate(start(6), step(6), best(6))
      call svensson_grid_start(rate, maturity, lower_tau, upper_tau, objective_name, start)
      start(5:6) = log(start(5:6))
      step(1:4) = max(0.05_dp * max(1.0_dp, abs(start(1:4))), 1.0e-4_dp)
      step(5:6) = 0.15_dp
      maxit = 3000
      if (present(max_iterations)) maxit = max(1, max_iterations)
      tol = 1.0e-10_dp
      if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
      call optimize_term_model(model_svensson, objective_name, rate, maturity, start, step, &
         log(lower_tau), log(upper_tau), best, fbest, fit%converged, &
         fit%iterations, evals, maxit, tol)
      fit%evaluations = evals
      fit%parameters(1:4) = best(1:4)
      fit%parameters(5:6) = exp(best(5:6))
      fit%objective_value = fbest
      fit%fitted = svensson_curve(maturity, fit%parameters)
      call finish_fit(rate, fit)
   end subroutine fit_svensson

   subroutine initialize_fit(fit, model_name, objective_name, n, npar)
      type(term_structure_fit), intent(out) :: fit
      character(len=*), intent(in) :: model_name, objective_name
      integer, intent(in) :: n, npar
      fit%model = model_name
      fit%objective = objective_name
      allocate(fit%parameters(npar), fit%fitted(n), fit%residuals(n))
      fit%parameters = 0.0_dp
      fit%fitted = 0.0_dp
      fit%residuals = 0.0_dp
   end subroutine initialize_fit

   pure logical function valid_data(rate, maturity, min_n)
      real(dp), intent(in) :: rate(:), maturity(:)
      integer, intent(in) :: min_n
      valid_data = size(rate) == size(maturity) .and. size(rate) >= min_n .and. &
         all(ieee_is_finite(rate)) .and. all(ieee_is_finite(maturity)) .and. &
         all(maturity >= 0.0_dp) .and. any(maturity > 0.0_dp)
   end function valid_data

   subroutine choose_tau_bounds(maturity, lower_in, upper_in, lower_tau, upper_tau)
      real(dp), intent(in) :: maturity(:)
      real(dp), intent(in), optional :: lower_in, upper_in
      real(dp), intent(out) :: lower_tau, upper_tau
      real(dp) :: min_positive, max_maturity
      min_positive = minval(maturity, mask=maturity > 0.0_dp)
      max_maturity = maxval(maturity)
      lower_tau = max(1.0e-4_dp, 0.05_dp * min_positive)
      upper_tau = max(10.0_dp * max_maturity, 100.0_dp * lower_tau)
      if (present(lower_in)) lower_tau = max(lower_in, 1.0e-8_dp)
      if (present(upper_in)) upper_tau = max(upper_in, lower_tau * (1.0_dp + 1.0e-6_dp))
   end subroutine choose_tau_bounds

   subroutine ns_design(maturity, tau, x)
      real(dp), intent(in) :: maturity(:), tau
      real(dp), intent(out) :: x(:, :)
      real(dp) :: z(size(maturity))
      z = maturity / tau
      x(:, 1) = 1.0_dp
      x(:, 2) = ns_loading1(z)
      x(:, 3) = ns_loading2(z)
   end subroutine ns_design

   subroutine svensson_design(maturity, tau1, tau2, x)
      real(dp), intent(in) :: maturity(:), tau1, tau2
      real(dp), intent(out) :: x(:, :)
      real(dp) :: z1(size(maturity)), z2(size(maturity))
      z1 = maturity / tau1
      z2 = maturity / tau2
      x(:, 1) = 1.0_dp
      x(:, 2) = exp(-z1)
      x(:, 3) = z1 * exp(-z1)
      x(:, 4) = z2 * exp(-z2)
   end subroutine svensson_design

   subroutine ns_grid_start(rate, maturity, lower_tau, upper_tau, start)
      real(dp), intent(in) :: rate(:), maturity(:), lower_tau, upper_tau
      real(dp), intent(out) :: start(4)
      real(dp), allocatable :: x(:, :), beta(:), fitted(:)
      real(dp) :: tau, value, best_value, log_lower, log_upper
      integer :: i, rank, info, ngrid
      allocate(x(size(rate), 3), beta(3), fitted(size(rate)))
      best_value = huge(1.0_dp)
      start = 0.0_dp
      ngrid = max(64, size(rate))
      log_lower = log(lower_tau)
      log_upper = log(upper_tau)
      do i = 1, ngrid
         tau = exp(log_lower + real(i - 1, dp) * (log_upper - log_lower) / real(ngrid - 1, dp))
         call ns_design(maturity, tau, x)
         call least_squares(x, rate, beta, rank, info)
         if (info /= 0 .or. rank < 3) cycle
         fitted = matmul(x, beta)
         value = sum((rate - fitted)**2)
         if (value < best_value) then
            best_value = value
            start(1:3) = beta
            start(4) = tau
         end if
      end do
   end subroutine ns_grid_start

   subroutine svensson_grid_start(rate, maturity, lower_tau, upper_tau, objective, start)
      real(dp), intent(in) :: rate(:), maturity(:), lower_tau, upper_tau
      character(len=*), intent(in) :: objective
      real(dp), intent(out) :: start(6)
      real(dp), allocatable :: x(:, :), beta(:), fitted(:), grid(:)
      real(dp) :: value, best_value
      integer :: i, j, rank, info, ngrid
      allocate(x(size(rate), 4), beta(4), fitted(size(rate)))
      ngrid = min(48, max(16, size(rate)))
      allocate(grid(ngrid))
      do i = 1, ngrid
         grid(i) = exp(log(lower_tau) + real(i - 1, dp) * &
            (log(upper_tau) - log(lower_tau)) / real(ngrid - 1, dp))
      end do
      best_value = huge(1.0_dp)
      start = 0.0_dp
      do i = 1, ngrid
         do j = 1, ngrid
            if (abs(log(grid(i) / grid(j))) < 1.0e-5_dp) cycle
            call svensson_design(maturity, grid(i), grid(j), x)
            call least_squares(x, rate, beta, rank, info)
            if (info /= 0 .or. rank < 4) cycle
            fitted = matmul(x, beta)
            value = loss_value(rate - fitted, objective)
            if (value < best_value) then
               best_value = value
               start(1:4) = beta
               start(5) = grid(i)
               start(6) = grid(j)
            end if
         end do
      end do
   end subroutine svensson_grid_start

   subroutine optimize_term_model(model, objective, rate, maturity, start, step, &
      log_tau_lower, log_tau_upper, best, fbest, converged, iterations, evaluations, &
      max_iterations, tolerance)
      integer, intent(in) :: model, max_iterations
      character(len=*), intent(in) :: objective
      real(dp), intent(in) :: rate(:), maturity(:), start(:), step(:)
      real(dp), intent(in) :: log_tau_lower, log_tau_upper, tolerance
      real(dp), intent(out) :: best(:), fbest
      logical, intent(out) :: converged
      integer, intent(out) :: iterations, evaluations
      real(dp), allocatable :: simplex(:, :), values(:), centroid(:), xr(:), xe(:), xc(:)
      real(dp) :: fr, fe, fc, diameter
      integer :: n, j

      n = size(start)
      allocate(simplex(n, n + 1), values(n + 1), centroid(n), xr(n), xe(n), xc(n))
      simplex(:, 1) = start
      do j = 1, n
         simplex(:, j + 1) = start
         simplex(j, j + 1) = simplex(j, j + 1) + step(j)
      end do
      evaluations = 0
      do j = 1, n + 1
         values(j) = model_objective(model, objective, simplex(:, j), rate, maturity, &
            log_tau_lower, log_tau_upper)
         evaluations = evaluations + 1
      end do
      converged = .false.
      do iterations = 1, max_iterations
         call sort_simplex(simplex, values)
         diameter = maxval(abs(simplex(:, 2:n + 1) - spread(simplex(:, 1), 2, n)))
         if (diameter <= tolerance * (1.0_dp + maxval(abs(simplex(:, 1)))) .and. &
             maxval(abs(values(2:n + 1) - values(1))) <= &
             tolerance * (1.0_dp + abs(values(1)))) then
            converged = .true.
            exit
         end if
         centroid = sum(simplex(:, 1:n), dim=2) / real(n, dp)
         xr = centroid + (centroid - simplex(:, n + 1))
         fr = model_objective(model, objective, xr, rate, maturity, log_tau_lower, log_tau_upper)
         evaluations = evaluations + 1
         if (fr < values(1)) then
            xe = centroid + 2.0_dp * (xr - centroid)
            fe = model_objective(model, objective, xe, rate, maturity, log_tau_lower, log_tau_upper)
            evaluations = evaluations + 1
            if (fe < fr) then
               simplex(:, n + 1) = xe
               values(n + 1) = fe
            else
               simplex(:, n + 1) = xr
               values(n + 1) = fr
            end if
         else if (fr < values(n)) then
            simplex(:, n + 1) = xr
            values(n + 1) = fr
         else
            if (fr < values(n + 1)) then
               xc = centroid + 0.5_dp * (xr - centroid)
            else
               xc = centroid + 0.5_dp * (simplex(:, n + 1) - centroid)
            end if
            fc = model_objective(model, objective, xc, rate, maturity, log_tau_lower, log_tau_upper)
            evaluations = evaluations + 1
            if (fc < min(fr, values(n + 1))) then
               simplex(:, n + 1) = xc
               values(n + 1) = fc
            else
               do j = 2, n + 1
                  simplex(:, j) = simplex(:, 1) + 0.5_dp * (simplex(:, j) - simplex(:, 1))
                  values(j) = model_objective(model, objective, simplex(:, j), rate, maturity, &
                     log_tau_lower, log_tau_upper)
                  evaluations = evaluations + 1
               end do
            end if
         end if
      end do
      call sort_simplex(simplex, values)
      best = simplex(:, 1)
      fbest = values(1)
      if (iterations > max_iterations) iterations = max_iterations
   end subroutine optimize_term_model

   function model_objective(model, objective, z, rate, maturity, log_tau_lower, &
      log_tau_upper) result(value)
      integer, intent(in) :: model
      character(len=*), intent(in) :: objective
      real(dp), intent(in) :: z(:), rate(:), maturity(:), log_tau_lower, log_tau_upper
      real(dp) :: value
      real(dp), allocatable :: parameters(:), fitted(:)
      integer :: first_tau
      first_tau = merge(4, 5, model == model_ns)
      if (any(z(first_tau:) < log_tau_lower) .or. any(z(first_tau:) > log_tau_upper)) then
         value = huge(1.0_dp) / 100.0_dp
         return
      end if
      allocate(parameters(size(z)), fitted(size(rate)))
      parameters = z
      parameters(first_tau:) = exp(z(first_tau:))
      select case (model)
      case (model_ns)
         fitted = nelson_siegel_curve(maturity, parameters)
      case (model_svensson)
         fitted = svensson_curve(maturity, parameters)
      case default
         value = huge(1.0_dp)
         return
      end select
      if (.not. all(ieee_is_finite(fitted))) then
         value = huge(1.0_dp) / 100.0_dp
      else
         value = loss_value(rate - fitted, objective)
      end if
   end function model_objective

   pure function loss_value(residuals, objective) result(value)
      real(dp), intent(in) :: residuals(:)
      character(len=*), intent(in) :: objective
      real(dp) :: value
      if (trim(objective) == "l1") then
         value = sum(abs(residuals))
      else
         value = sum(residuals**2)
      end if
   end function loss_value

   subroutine finish_fit(rate, fit)
      real(dp), intent(in) :: rate(:)
      type(term_structure_fit), intent(inout) :: fit
      fit%residuals = rate - fit%fitted
      fit%sse = sum(fit%residuals**2)
      fit%mae = sum(abs(fit%residuals)) / real(size(rate), dp)
      fit%rmse = sqrt(fit%sse / real(size(rate), dp))
      fit%status = merge(0, 1, fit%converged)
   end subroutine finish_fit

   subroutine sort_simplex(points, values)
      real(dp), intent(inout) :: points(:, :), values(:)
      real(dp), allocatable :: point_tmp(:)
      real(dp) :: value_tmp
      integer :: i, j, k
      allocate(point_tmp(size(points, 1)))
      do i = 1, size(values) - 1
         k = i
         do j = i + 1, size(values)
            if (values(j) < values(k)) k = j
         end do
         if (k /= i) then
            value_tmp = values(i)
            values(i) = values(k)
            values(k) = value_tmp
            point_tmp = points(:, i)
            points(:, i) = points(:, k)
            points(:, k) = point_tmp
         end if
      end do
   end subroutine sort_simplex

   pure function lower_string(text) result(lower)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i, code
      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + 32)
      end do
   end function lower_string
end module fbonds_term_structure
