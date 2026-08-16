! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_timeseries
   use vgam_kinds, only : dp, pi
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   type, public :: ar1_result_t
      real(dp) :: drift = 0.0_dp
      real(dp) :: innovation_sd = 1.0_dp
      real(dp) :: rho = 0.0_dp
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
      logical :: conditional = .false.
   contains
      procedure :: forecast => forecast_ar1
   end type ar1_result_t

   public :: dar1_log, dar1_density, fit_ar1

contains

   subroutine dar1_log(x, drift, var_error, rho, log_density, conditional)
      real(dp), intent(in) :: x(:), drift, var_error, rho
      real(dp), allocatable, intent(out) :: log_density(:)
      logical, intent(in), optional :: conditional
      logical :: cond
      real(dp) :: mean0, var0
      integer :: i, n

      n = size(x)
      allocate(log_density(n))
      cond = .false.
      if (present(conditional)) cond = conditional
      if (n == 0 .or. var_error <= 0.0_dp .or. abs(rho) >= 1.0_dp) then
         log_density = -huge(1.0_dp)
         return
      end if
      mean0 = drift/(1.0_dp - rho)
      var0 = var_error/(1.0_dp - rho*rho)
      log_density(1) = normal_logpdf(x(1), mean0, var0)
      if (cond) log_density(1) = 0.0_dp
      do i = 2, n
         log_density(i) = normal_logpdf(x(i), drift + rho*x(i - 1), var_error)
      end do
   end subroutine dar1_log

   subroutine dar1_density(x, drift, var_error, rho, density, conditional)
      real(dp), intent(in) :: x(:), drift, var_error, rho
      real(dp), allocatable, intent(out) :: density(:)
      logical, intent(in), optional :: conditional
      call dar1_log(x, drift, var_error, rho, density, conditional)
      density = exp(density)
   end subroutine dar1_density

   subroutine fit_ar1(x, result, conditional, max_iter, tol)
      real(dp), intent(in) :: x(:)
      type(ar1_result_t), intent(out) :: result
      logical, intent(in), optional :: conditional
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: par(:), hess(:, :), cov(:, :), jac(:, :), covn(:, :)
      real(dp) :: fval, tolerance, rho0, v0, den
      integer :: n, stat, stat2, niter, i

      n = size(x)
      if (n < 3) then
         result%status = 1
         return
      end if
      result%conditional = .false.
      if (present(conditional)) result%conditional = conditional
      den = sum((x(1:n - 1) - sum(x(1:n - 1))/real(n - 1, dp))**2)
      if (den > tiny(1.0_dp)) then
         rho0 = sum((x(1:n - 1) - sum(x(1:n - 1))/real(n - 1, dp))* &
                    (x(2:n) - sum(x(2:n))/real(n - 1, dp)))/den
      else
         rho0 = 0.0_dp
      end if
      rho0 = min(0.95_dp, max(-0.95_dp, rho0))
      v0 = max(sum((x(2:n) - rho0*x(1:n - 1))**2)/real(n - 1, dp), 1.0e-8_dp)
      allocate(par(3))
      par(1) = (1.0_dp - rho0)*sum(x)/real(n, dp)
      par(2) = log(sqrt(v0))
      par(3) = atanh(rho0)
      niter = 250
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-8_dp
      if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, niter, tolerance)
      result%status = stat
      result%converged = stat == 0
      result%drift = par(1)
      result%innovation_sd = exp(min(par(2), 50.0_dp))
      result%rho = tanh(par(3))
      allocate(result%fitted(n), result%residuals(n))
      result%fitted(1) = result%drift/(1.0_dp - result%rho)
      do i = 2, n
         result%fitted(i) = result%drift + result%rho*x(i - 1)
      end do
      result%residuals = x - result%fitted
      if (result%conditional) result%residuals(1) = 0.0_dp
      result%loglik = -fval
      result%aic = 2.0_dp*fval + 6.0_dp
      allocate(hess(3, 3))
      call numerical_hessian(objective, par, hess)
      call invert_matrix(hess, cov, stat2)
      if (stat2 == 0) then
         allocate(jac(3, 3), covn(3, 3))
         jac = 0.0_dp
         jac(1, 1) = 1.0_dp
         jac(2, 2) = result%innovation_sd
         jac(3, 3) = 1.0_dp - result%rho*result%rho
         covn = matmul(jac, matmul(cov, transpose(jac)))
         result%covariance = covn
      else
         allocate(result%covariance(3, 3))
         result%covariance = 0.0_dp
      end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp), allocatable :: lden(:)
         real(dp) :: sd, rr
         sd = exp(min(theta(2), 50.0_dp))
         rr = tanh(theta(3))
         call dar1_log(x, theta(1), sd*sd, rr, lden, result%conditional)
         nll = -sum(lden)
      end function objective
   end subroutine fit_ar1

   subroutine forecast_ar1(self, last_value, steps, mean, variance)
      class(ar1_result_t), intent(in) :: self
      real(dp), intent(in) :: last_value
      integer, intent(in) :: steps
      real(dp), allocatable, intent(out) :: mean(:), variance(:)
      real(dp) :: level, powrho, acc
      integer :: h
      if (steps < 1) then
         allocate(mean(0), variance(0))
         return
      end if
      allocate(mean(steps), variance(steps))
      level = self%drift/(1.0_dp - self%rho)
      acc = 0.0_dp
      do h = 1, steps
         powrho = self%rho**h
         mean(h) = level + powrho*(last_value - level)
         acc = acc + self%rho**(2*(h - 1))
         variance(h) = self%innovation_sd**2*acc
      end do
   end subroutine forecast_ar1

   elemental real(dp) function normal_logpdf(x, mean, variance) result(v)
      real(dp), intent(in) :: x, mean, variance
      if (variance <= 0.0_dp) then
         v = -huge(1.0_dp)
      else
         v = -0.5_dp*(log(2.0_dp*pi*variance) + (x - mean)**2/variance)
      end if
   end function normal_logpdf

end module vgam_timeseries
