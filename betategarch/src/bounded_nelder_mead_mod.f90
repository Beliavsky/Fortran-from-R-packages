! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2011-2025 Genaro Sucarrat
! Copyright (C) 2026 contributors to the Modern Fortran translation
!
! This file is part of betategarch-modern-fortran.
! It is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License version 2 only.

module bounded_nelder_mead_mod
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use betategarch_kinds, only : dp
  implicit none
  private

  public :: bounded_nelder_mead
  public :: parameter_to_unconstrained, unconstrained_to_parameter

  abstract interface
    function objective_function(x, context) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      class(*), intent(in) :: context
      real(dp) :: value
    end function objective_function
  end interface

contains

  pure subroutine unconstrained_to_parameter(z, lower, upper, x)
    real(dp), intent(in) :: z(:), lower(:), upper(:)
    real(dp), intent(out) :: x(:)

    real(dp) :: ez, span
    logical :: has_lower, has_upper
    integer :: i

    do i = 1, size(z)
      has_lower = lower(i) > -sqrt(huge(1.0_dp))
      has_upper = upper(i) < sqrt(huge(1.0_dp))
      if (has_lower .and. has_upper) then
        span = upper(i) - lower(i)
        if (z(i) >= 0.0_dp) then
          ez = exp(-min(z(i), 700.0_dp))
          x(i) = lower(i) + span/(1.0_dp + ez)
        else
          ez = exp(max(z(i), -700.0_dp))
          x(i) = lower(i) + span*ez/(1.0_dp + ez)
        end if
      else if (has_lower) then
        x(i) = lower(i) + exp(min(z(i), 700.0_dp))
      else if (has_upper) then
        x(i) = upper(i) - exp(min(z(i), 700.0_dp))
      else
        x(i) = z(i)
      end if
    end do
  end subroutine unconstrained_to_parameter

  pure subroutine parameter_to_unconstrained(x, lower, upper, z)
    real(dp), intent(in) :: x(:), lower(:), upper(:)
    real(dp), intent(out) :: z(:)

    real(dp) :: ratio, delta, eps_bound
    logical :: has_lower, has_upper
    integer :: i

    eps_bound = sqrt(epsilon(1.0_dp))
    do i = 1, size(x)
      has_lower = lower(i) > -sqrt(huge(1.0_dp))
      has_upper = upper(i) < sqrt(huge(1.0_dp))
      if (has_lower .and. has_upper) then
        ratio = (x(i) - lower(i))/(upper(i) - lower(i))
        ratio = min(1.0_dp - eps_bound, max(eps_bound, ratio))
        z(i) = log(ratio/(1.0_dp - ratio))
      else if (has_lower) then
        delta = max(x(i) - lower(i), eps_bound)
        z(i) = log(delta)
      else if (has_upper) then
        delta = max(upper(i) - x(i), eps_bound)
        z(i) = log(delta)
      else
        z(i) = x(i)
      end if
    end do
  end subroutine parameter_to_unconstrained

  subroutine bounded_nelder_mead(fun, context, initial, lower, upper, solution, fval, convergence, &
      iterations, evaluations, max_iterations, tolerance, initial_step)
    procedure(objective_function) :: fun
    class(*), intent(in) :: context
    real(dp), intent(in) :: initial(:), lower(:), upper(:)
    real(dp), intent(out) :: solution(:)
    real(dp), intent(out) :: fval
    integer, intent(out) :: convergence, iterations, evaluations
    integer, intent(in), optional :: max_iterations
    real(dp), intent(in), optional :: tolerance, initial_step

    real(dp), allocatable :: simplex(:, :), values(:), centroid(:), xr(:), xe(:), xc(:)
    real(dp), allocatable :: z0(:), xwork(:)
    real(dp) :: alpha, gamma, rho, sigma, tol, step_scale
    real(dp) :: fr, fe, fc, fspread, xspread
    integer :: n, max_iter, i, iter

    n = size(initial)
    if (size(lower) /= n .or. size(upper) /= n .or. size(solution) /= n) then
      error stop "bounded_nelder_mead: inconsistent dimensions"
    end if

    max_iter = 2000
    if (present(max_iterations)) max_iter = max_iterations
    tol = 1.0e-7_dp
    if (present(tolerance)) tol = tolerance
    step_scale = 0.10_dp
    if (present(initial_step)) step_scale = initial_step

    alpha = 1.0_dp
    gamma = 2.0_dp
    rho = 0.5_dp
    sigma = 0.5_dp

    allocate(simplex(n, n+1), values(n+1), centroid(n), xr(n), xe(n), xc(n), z0(n), xwork(n))
    call parameter_to_unconstrained(initial, lower, upper, z0)
    simplex(:, 1) = z0
    do i = 1, n
      simplex(:, i+1) = z0
      simplex(i, i+1) = simplex(i, i+1) + step_scale*max(1.0_dp, abs(z0(i)))
    end do

    evaluations = 0
    do i = 1, n + 1
      values(i) = evaluate_z(simplex(:, i))
    end do
    call sort_simplex(simplex, values)

    convergence = 1
    iterations = 0
    do iter = 1, max_iter
      iterations = iter
      fspread = maxval(abs(values - values(1)))
      xspread = maxval(abs(simplex - spread(simplex(:, 1), dim=2, ncopies=n+1)))
      if (fspread <= tol*(1.0_dp + abs(values(1))) .and. &
          xspread <= sqrt(tol)*(1.0_dp + maxval(abs(simplex(:, 1))))) then
        convergence = 0
        exit
      end if

      centroid = sum(simplex(:, 1:n), dim=2)/real(n, dp)
      xr = centroid + alpha*(centroid - simplex(:, n+1))
      fr = evaluate_z(xr)

      if (fr < values(1)) then
        xe = centroid + gamma*(xr - centroid)
        fe = evaluate_z(xe)
        if (fe < fr) then
          simplex(:, n+1) = xe
          values(n+1) = fe
        else
          simplex(:, n+1) = xr
          values(n+1) = fr
        end if
      else if (fr < values(n)) then
        simplex(:, n+1) = xr
        values(n+1) = fr
      else
        if (fr < values(n+1)) then
          xc = centroid + rho*(xr - centroid)
        else
          xc = centroid + rho*(simplex(:, n+1) - centroid)
        end if
        fc = evaluate_z(xc)
        if (fc < min(fr, values(n+1))) then
          simplex(:, n+1) = xc
          values(n+1) = fc
        else
          do i = 2, n + 1
            simplex(:, i) = simplex(:, 1) + sigma*(simplex(:, i) - simplex(:, 1))
            values(i) = evaluate_z(simplex(:, i))
          end do
        end if
      end if
      call sort_simplex(simplex, values)
    end do

    call unconstrained_to_parameter(simplex(:, 1), lower, upper, solution)
    fval = values(1)

  contains

    function evaluate_z(z) result(value)
      real(dp), intent(in) :: z(:)
      real(dp) :: value

      call unconstrained_to_parameter(z, lower, upper, xwork)
      value = fun(xwork, context)
      if (.not. ieee_is_finite(value)) value = huge(1.0_dp)/100.0_dp
      evaluations = evaluations + 1
    end function evaluate_z

  end subroutine bounded_nelder_mead

  subroutine sort_simplex(simplex, values)
    real(dp), intent(inout) :: simplex(:, :), values(:)

    real(dp), allocatable :: column(:)
    real(dp) :: value_tmp
    integer :: i, j, best

    allocate(column(size(simplex, 1)))
    do i = 1, size(values) - 1
      best = i
      do j = i + 1, size(values)
        if (values(j) < values(best)) best = j
      end do
      if (best /= i) then
        value_tmp = values(i)
        values(i) = values(best)
        values(best) = value_tmp
        column = simplex(:, i)
        simplex(:, i) = simplex(:, best)
        simplex(:, best) = column
      end if
    end do
  end subroutine sort_simplex

end module bounded_nelder_mead_mod
