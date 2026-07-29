! SPDX-License-Identifier: GPL-3.0-only
!
! Copyright (C) 2017 Bernhard Pfaff
! Modern Fortran translation of the computational algorithms in mcrp.
module mcrp_optimizer
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use mcrp_kinds, only : dp, mcrp_success, mcrp_invalid_shape, &
      mcrp_invalid_argument, mcrp_numerical_failure, mcrp_max_iterations
   use mcrp_moments, only : m2, port_risk_contrib, port_skew_contrib, &
      port_kurt_contrib
   implicit none
   private

   type, public :: mcrp_result
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: raw_parameters(:)
      real(dp), allocatable :: variance_contributions(:)
      real(dp), allocatable :: skewness_contributions(:)
      real(dp), allocatable :: kurtosis_contributions(:)
      real(dp) :: objective = huge(1.0_dp)
      integer :: iterations = 0
      integer :: evaluations = 0
      integer :: status = mcrp_invalid_argument
      logical :: converged = .false.
      character(len=160) :: message = ''
   end type mcrp_result

   public :: mcrp
   public :: mcrp_objective_value


contains

   subroutine mcrp(start, returns, result, lambda, active, lower, upper, &
      max_iterations, tolerance, initial_step)
      real(dp), intent(in) :: start(:), returns(:, :)
      type(mcrp_result), intent(out) :: result
      real(dp), intent(in), optional :: lambda(3)
      logical, intent(in), optional :: active(3)
      real(dp), intent(in), optional :: lower(:), upper(:)
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance, initial_step

      real(dp), allocatable :: lo(:), hi(:), rc(:, :), covariance(:, :)
      real(dp), allocatable :: best(:)
      real(dp) :: lam(3), tol, step, scale
      logical :: enabled(3)
      integer :: nasset, nobs, maxit, opt_status, iterations, evaluations
      integer :: i

      nobs = size(returns, 1)
      nasset = size(returns, 2)
      if (nobs < 2 .or. nasset < 1 .or. size(start) /= nasset) then
         call fail_result(result, mcrp_invalid_shape, &
            'returns and start have incompatible dimensions')
         return
      end if
      if (.not. all(ieee_is_finite(returns)) .or. &
          .not. all(ieee_is_finite(start))) then
         call fail_result(result, mcrp_invalid_argument, &
            'inputs must be finite')
         return
      end if

      lam = 1.0_dp
      if (present(lambda)) lam = lambda
      enabled = .not. ieee_is_nan(lam)
      if (present(active)) enabled = enabled .and. active
      do i = 1, 3
         if (.not. enabled(i)) lam(i) = 0.0_dp
      end do
      if (.not. any(enabled)) then
         call fail_result(result, mcrp_invalid_argument, &
            'at least one criterion must be active')
         return
      end if

      allocate(lo(nasset), hi(nasset))
      lo = -huge(1.0_dp) / 100.0_dp
      hi = huge(1.0_dp) / 100.0_dp
      if (present(lower)) then
         if (size(lower) /= nasset) then
            call fail_result(result, mcrp_invalid_shape, &
               'lower bound has the wrong size')
            return
         end if
         lo = lower
      end if
      if (present(upper)) then
         if (size(upper) /= nasset) then
            call fail_result(result, mcrp_invalid_shape, &
               'upper bound has the wrong size')
            return
         end if
         hi = upper
      end if
      if (any(lo > hi) .or. any(start < lo) .or. any(start > hi)) then
         call fail_result(result, mcrp_invalid_argument, &
            'invalid bounds or starting values outside bounds')
         return
      end if

      maxit = 4000
      if (present(max_iterations)) maxit = max_iterations
      tol = 1.0e-10_dp
      if (present(tolerance)) tol = tolerance
      step = 0.10_dp
      if (present(initial_step)) step = initial_step
      if (maxit < 1 .or. tol <= 0.0_dp .or. step <= 0.0_dp) then
         call fail_result(result, mcrp_invalid_argument, &
            'invalid optimizer controls')
         return
      end if

      rc = center_data(returns)
      covariance = m2(returns)

      call bounded_nelder_mead(start, lo, hi, maxit, tol, step, covariance, &
         rc, lam, enabled, best, result%objective, iterations, evaluations, &
         opt_status)

      scale = sum(abs(best))
      if (.not. ieee_is_finite(scale) .or. scale <= tiny(1.0_dp)) then
         call fail_result(result, mcrp_numerical_failure, &
            'optimizer returned a zero or nonfinite parameter vector')
         return
      end if

      result%raw_parameters = best
      result%weights = best / scale
      result%iterations = iterations
      result%evaluations = evaluations
      result%status = opt_status
      result%converged = opt_status == mcrp_success
      if (result%converged) then
         result%message = 'converged'
      else
         result%message = 'maximum iterations reached; best point returned'
      end if

      if (enabled(1)) then
         result%variance_contributions = &
            port_risk_contrib(returns, result%weights, .true.)
      else
         allocate(result%variance_contributions(0))
      end if
      if (enabled(2)) then
         result%skewness_contributions = &
            port_skew_contrib(returns, result%weights, .true.)
      else
         allocate(result%skewness_contributions(0))
      end if
      if (enabled(3)) then
         result%kurtosis_contributions = &
            port_kurt_contrib(returns, result%weights, .true.)
      else
         allocate(result%kurtosis_contributions(0))
      end if

   end subroutine mcrp

   function mcrp_objective_value(returns, x, lambda, active, status) result(value)
      real(dp), intent(in) :: returns(:, :), x(:)
      real(dp), intent(in), optional :: lambda(3)
      logical, intent(in), optional :: active(3)
      integer, intent(out), optional :: status
      real(dp) :: value
      real(dp), allocatable :: covariance(:, :), rc(:, :)
      real(dp) :: lam(3)
      logical :: enabled(3)
      integer :: i

      if (size(returns, 1) < 2 .or. size(returns, 2) /= size(x) .or. &
          size(x) < 1) then
         value = huge(1.0_dp) / 1000.0_dp
         if (present(status)) status = mcrp_invalid_shape
         return
      end if
      lam = 1.0_dp
      if (present(lambda)) lam = lambda
      enabled = .not. ieee_is_nan(lam)
      if (present(active)) enabled = enabled .and. active
      do i = 1, 3
         if (.not. enabled(i)) lam(i) = 0.0_dp
      end do
      if (.not. any(enabled)) then
         value = huge(1.0_dp) / 1000.0_dp
         if (present(status)) status = mcrp_invalid_argument
         return
      end if
      covariance = m2(returns)
      rc = center_data(returns)
      value = objective_from_workspace(x, covariance, rc, lam, enabled)
      if (present(status)) then
         if (value >= huge(1.0_dp) / 1.0e4_dp) then
            status = mcrp_numerical_failure
         else
            status = mcrp_success
         end if
      end if
   end function mcrp_objective_value

   function objective_from_workspace(x, covariance, rc, lambda, active) &
      result(value)
      real(dp), intent(in) :: x(:), covariance(:, :), rc(:, :), lambda(3)
      logical, intent(in) :: active(3)
      real(dp) :: value
      real(dp), allocatable :: covx(:), s(:), deriv(:), pct(:)
      real(dp) :: variance, third, fourth, skewness, kurtosis
      real(dp), parameter :: penalty = huge(1.0_dp) / 1000.0_dp
      integer :: nobs

      if (.not. all(ieee_is_finite(x)) .or. &
          sum(abs(x)) <= sqrt(tiny(1.0_dp))) then
         value = penalty
         return
      end if

      nobs = size(rc, 1)
      covx = matmul(covariance, x)
      variance = dot_product(x, covx)
      if (.not. ieee_is_finite(variance) .or. variance <= tiny(1.0_dp)) then
         value = penalty
         return
      end if

      value = 0.0_dp
      if (active(1)) then
         pct = x * covx / (2.0_dp * variance)
         value = value + lambda(1) * sample_variance(pct)
      end if

      if (active(2) .or. active(3)) s = matmul(rc, x)

      if (active(2)) then
         third = sum(s**3) / real(nobs, dp)
         skewness = third / variance**1.5_dp
         if (.not. ieee_is_finite(skewness) .or. &
             abs(skewness) <= sqrt(tiny(1.0_dp))) then
            value = penalty
            return
         end if
         deriv = 3.0_dp * matmul(transpose(rc), s**2) / real(nobs, dp)
         deriv = (variance**1.5_dp * deriv - &
            2.0_dp * third * sqrt(variance) * covx) / variance**3
         pct = x * deriv / skewness
         value = value + lambda(2) * sample_variance(pct)
      end if

      if (active(3)) then
         fourth = sum(s**4) / real(nobs, dp)
         kurtosis = fourth / variance**2
         if (.not. ieee_is_finite(kurtosis) .or. &
             abs(kurtosis) <= sqrt(tiny(1.0_dp))) then
            value = penalty
            return
         end if
         deriv = 4.0_dp * matmul(transpose(rc), s**3) / real(nobs, dp)
         deriv = (variance * deriv - 2.0_dp * fourth * covx) / &
            (2.0_dp * variance**3)
         pct = x * deriv / kurtosis
         value = value + lambda(3) * sample_variance(pct)
      end if

      if (.not. ieee_is_finite(value)) value = penalty
   end function objective_from_workspace

   subroutine bounded_nelder_mead(start, lower, upper, max_iterations, &
      tolerance, initial_step, covariance, rc, lambda, active, best, fbest, &
      iterations, evaluations, status)
      real(dp), intent(in) :: start(:), lower(:), upper(:)
      integer, intent(in) :: max_iterations
      real(dp), intent(in) :: tolerance, initial_step
      real(dp), intent(in) :: covariance(:, :), rc(:, :), lambda(3)
      logical, intent(in) :: active(3)
      real(dp), allocatable, intent(out) :: best(:)
      real(dp), intent(out) :: fbest
      integer, intent(out) :: iterations, evaluations, status

      real(dp), allocatable :: simplex(:, :), f(:), centroid(:), xr(:), xe(:)
      real(dp), allocatable :: xc(:)
      real(dp) :: fr, fe, fc, delta, fspread, diameter
      real(dp), parameter :: alpha = 1.0_dp, gamma = 2.0_dp
      real(dp), parameter :: rho = 0.5_dp, sigma = 0.5_dp
      integer :: n, i, j

      n = size(start)
      allocate(simplex(n, n + 1), f(n + 1), centroid(n), xr(n), xe(n), xc(n))
      simplex(:, 1) = project_bounds(start, lower, upper)
      do j = 1, n
         simplex(:, j + 1) = simplex(:, 1)
         delta = initial_step * max(1.0_dp, abs(start(j)))
         simplex(j, j + 1) = min(upper(j), simplex(j, 1) + delta)
         if (abs(simplex(j, j + 1) - simplex(j, 1)) <= tiny(1.0_dp)) then
            simplex(j, j + 1) = max(lower(j), simplex(j, 1) - delta)
         end if
         if (abs(simplex(j, j + 1) - simplex(j, 1)) <= tiny(1.0_dp) .and. &
             lower(j) < upper(j)) then
            simplex(j, j + 1) = 0.5_dp * (lower(j) + upper(j))
         end if
      end do

      do i = 1, n + 1
         f(i) = objective_from_workspace(simplex(:, i), covariance, rc, lambda, active)
      end do
      evaluations = n + 1
      status = mcrp_max_iterations

      do iterations = 1, max_iterations
         call sort_simplex(simplex, f)
         fspread = maxval(abs(f(2:) - f(1)))
         diameter = 0.0_dp
         do i = 2, n + 1
            diameter = max(diameter, vector_norm(simplex(:, i) - simplex(:, 1)))
         end do
         if (fspread <= tolerance * (1.0_dp + abs(f(1))) .and. &
             diameter <= sqrt(tolerance) * &
             (1.0_dp + vector_norm(simplex(:, 1)))) then
            status = mcrp_success
            exit
         end if

         centroid = sum(simplex(:, 1:n), dim=2) / real(n, dp)
         xr = project_bounds(centroid + alpha * &
            (centroid - simplex(:, n + 1)), lower, upper)
         fr = objective_from_workspace(xr, covariance, rc, lambda, active)
         evaluations = evaluations + 1

         if (fr < f(1)) then
            xe = project_bounds(centroid + gamma * (xr - centroid), lower, upper)
            fe = objective_from_workspace(xe, covariance, rc, lambda, active)
            evaluations = evaluations + 1
            if (fe < fr) then
               simplex(:, n + 1) = xe
               f(n + 1) = fe
            else
               simplex(:, n + 1) = xr
               f(n + 1) = fr
            end if
         else if (fr < f(n)) then
            simplex(:, n + 1) = xr
            f(n + 1) = fr
         else
            if (fr < f(n + 1)) then
               xc = project_bounds(centroid + rho * (xr - centroid), &
                  lower, upper)
            else
               xc = project_bounds(centroid - rho * (centroid - &
                  simplex(:, n + 1)), lower, upper)
            end if
            fc = objective_from_workspace(xc, covariance, rc, lambda, active)
            evaluations = evaluations + 1
            if (fc < min(fr, f(n + 1))) then
               simplex(:, n + 1) = xc
               f(n + 1) = fc
            else
               do i = 2, n + 1
                  simplex(:, i) = project_bounds(simplex(:, 1) + sigma * &
                     (simplex(:, i) - simplex(:, 1)), lower, upper)
                  f(i) = objective_from_workspace(simplex(:, i), covariance, rc, lambda, active)
               end do
               evaluations = evaluations + n
            end if
         end if
      end do

      call sort_simplex(simplex, f)
      best = simplex(:, 1)
      fbest = f(1)
      if (iterations > max_iterations) iterations = max_iterations
   end subroutine bounded_nelder_mead

   subroutine sort_simplex(simplex, f)
      real(dp), intent(inout) :: simplex(:, :), f(:)
      real(dp), allocatable :: tmp(:)
      real(dp) :: tf
      integer :: i, j, k

      allocate(tmp(size(simplex, 1)))
      do i = 1, size(f) - 1
         k = i
         do j = i + 1, size(f)
            if (f(j) < f(k)) k = j
         end do
         if (k /= i) then
            tf = f(i)
            f(i) = f(k)
            f(k) = tf
            tmp = simplex(:, i)
            simplex(:, i) = simplex(:, k)
            simplex(:, k) = tmp
         end if
      end do
   end subroutine sort_simplex

   function project_bounds(x, lower, upper) result(projected)
      real(dp), intent(in) :: x(:), lower(:), upper(:)
      real(dp), allocatable :: projected(:)
      projected = min(upper, max(lower, x))
   end function project_bounds

   function center_data(r) result(rc)
      real(dp), intent(in) :: r(:, :)
      real(dp), allocatable :: rc(:, :), means(:)
      integer :: j

      means = sum(r, dim=1) / real(size(r, 1), dp)
      allocate(rc(size(r, 1), size(r, 2)))
      do j = 1, size(r, 2)
         rc(:, j) = r(:, j) - means(j)
      end do
   end function center_data

   real(dp) function sample_variance(x)
      real(dp), intent(in) :: x(:)
      real(dp) :: mean_x

      if (size(x) <= 1) then
         sample_variance = 0.0_dp
      else
         mean_x = sum(x) / real(size(x), dp)
         sample_variance = sum((x - mean_x)**2) / real(size(x) - 1, dp)
      end if
   end function sample_variance

   real(dp) function vector_norm(x)
      real(dp), intent(in) :: x(:)
      vector_norm = sqrt(max(0.0_dp, dot_product(x, x)))
   end function vector_norm

   subroutine fail_result(result, status, message)
      type(mcrp_result), intent(out) :: result
      integer, intent(in) :: status
      character(len=*), intent(in) :: message
      result%status = status
      result%converged = .false.
      result%message = message
      allocate(result%weights(0), result%raw_parameters(0))
      allocate(result%variance_contributions(0))
      allocate(result%skewness_contributions(0))
      allocate(result%kurtosis_contributions(0))
   end subroutine fail_result

end module mcrp_optimizer
