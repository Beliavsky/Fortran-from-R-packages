! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_zero_altered_models
   use vgam_kinds, only : dp
   use vgam_zero_altered, only : dzapois_v, dzanbinom_v, dzageom_v, dzabinom_v
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   integer, parameter, public :: za_family_poisson = 1
   integer, parameter, public :: za_family_negbinomial = 2
   integer, parameter, public :: za_family_geometric = 3
   integer, parameter, public :: za_family_binomial = 4

   type, public :: zero_altered_count_result_t
      real(dp), allocatable :: parent_coefficients(:)
      real(dp), allocatable :: zero_coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_parent_parameter(:)
      real(dp), allocatable :: fitted_zero_probability(:)
      real(dp), allocatable :: fitted_mean(:)
      real(dp) :: size = huge(1.0_dp)
      integer :: trials = 0
      integer :: family = 0
      real(dp) :: loglik = -huge(1.0_dp), aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict => predict_zero_altered_count
   end type zero_altered_count_result_t

   public :: fit_zero_altered_poisson, fit_zero_altered_negative_binomial
   public :: fit_zero_altered_geometric, fit_zero_altered_binomial

contains

   subroutine fit_zero_altered_poisson(y, x_parent, x_zero, result, weights, max_iter, tol)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x_parent(:, :), x_zero(:, :)
      type(zero_altered_count_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      call fit_za_common(y, x_parent, x_zero, za_family_poisson, 0, result, weights, max_iter, tol)
   end subroutine fit_zero_altered_poisson

   subroutine fit_zero_altered_negative_binomial(y, x_parent, x_zero, result, weights, max_iter, tol)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x_parent(:, :), x_zero(:, :)
      type(zero_altered_count_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      call fit_za_common(y, x_parent, x_zero, za_family_negbinomial, 0, result, weights, max_iter, tol)
   end subroutine fit_zero_altered_negative_binomial

   subroutine fit_zero_altered_geometric(y, x_parent, x_zero, result, weights, max_iter, tol)
      integer, intent(in) :: y(:)
      real(dp), intent(in) :: x_parent(:, :), x_zero(:, :)
      type(zero_altered_count_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      call fit_za_common(y, x_parent, x_zero, za_family_geometric, 0, result, weights, max_iter, tol)
   end subroutine fit_zero_altered_geometric

   subroutine fit_zero_altered_binomial(y, trials, x_parent, x_zero, result, weights, max_iter, tol)
      integer, intent(in) :: y(:), trials
      real(dp), intent(in) :: x_parent(:, :), x_zero(:, :)
      type(zero_altered_count_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      call fit_za_common(y, x_parent, x_zero, za_family_binomial, trials, result, weights, max_iter, tol)
   end subroutine fit_zero_altered_binomial

   subroutine fit_za_common(y, xp, xz, family, trials, result, weights, max_iter, tol)
      integer, intent(in) :: y(:), family, trials
      real(dp), intent(in) :: xp(:, :), xz(:, :)
      type(zero_altered_count_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), h(:, :), cov(:, :)
      real(dp) :: fval, tolerance, mean_y, zfrac
      integer :: n, pp, pz, np, stat, stat2, niter, i

      n = size(y); pp = size(xp, 2); pz = size(xz, 2)
      if (n <= 0 .or. min(pp, pz) <= 0 .or. size(xp, 1) /= n .or. size(xz, 1) /= n .or. any(y < 0)) then
         result%status = 1; return
      end if
      if (family == za_family_binomial .and. (trials < 1 .or. any(y > trials))) then
         result%status = 2; return
      end if
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then; result%status = 3; return; end if
         w = weights
      end if
      np = pp + pz + merge(1, 0, family == za_family_negbinomial)
      allocate(par(np)); par = 0.0_dp
      mean_y = sum(w*real(y, dp))/max(sum(w), tiny(1.0_dp))
      zfrac = sum(w*merge(1.0_dp, 0.0_dp, y == 0))/max(sum(w), tiny(1.0_dp))
      if (is_intercept(xp)) then
         select case (family)
         case (za_family_poisson, za_family_negbinomial)
            par(1) = log(max(mean_y, 0.2_dp))
         case (za_family_geometric)
            par(1) = logit(min(0.95_dp, max(0.05_dp, 1.0_dp/max(1.0_dp, mean_y))))
         case (za_family_binomial)
            par(1) = logit(min(0.95_dp, max(0.05_dp, mean_y/real(trials, dp))))
         end select
      end if
      if (is_intercept(xz)) par(pp + 1) = logit(min(0.95_dp, max(0.05_dp, zfrac)))
      if (family == za_family_negbinomial) par(np) = log(2.0_dp)
      niter = 350; if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp; if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)

      result%status = stat; result%converged = stat == 0; result%family = family; result%trials = trials
      result%parent_coefficients = par(1:pp); result%zero_coefficients = par(pp + 1:pp + pz)
      if (family == za_family_negbinomial) result%size = exp(min(30.0_dp, par(np)))
      result%loglik = -fval; result%aic = 2.0_dp*fval + 2.0_dp*real(np, dp)
      allocate(result%fitted_parent_parameter(n), result%fitted_zero_probability(n), result%fitted_mean(n))
      do i = 1, n
         call fitted_values_row(result, xp(i, :), xz(i, :), result%fitted_parent_parameter(i), &
            result%fitted_zero_probability(i), result%fitted_mean(i))
      end do
      allocate(h(np, np)); call numerical_hessian(objective, par, h); call invert_matrix(h, cov, stat2)
      if (stat2 == 0) then; result%covariance = cov; else; allocate(result%covariance(0, 0)); end if
   contains
      real(dp) function objective(theta) result(nll)
         real(dp), intent(in) :: theta(:)
         real(dp) :: parent, pzero, sizev, d
         integer :: row
         nll = 0.0_dp; sizev = huge(1.0_dp)
         if (family == za_family_negbinomial) sizev = exp(min(30.0_dp, theta(np)))
         do row = 1, n
            if (family == za_family_poisson .or. family == za_family_negbinomial) then
               parent = exp(clamp_eta(dot_product(xp(row, :), theta(1:pp))))
            else
               parent = logistic(dot_product(xp(row, :), theta(1:pp)))
            end if
            pzero = logistic(dot_product(xz(row, :), theta(pp + 1:pp + pz)))
            select case (family)
            case (za_family_poisson)
               d = dzapois_v(y(row), parent, pzero)
            case (za_family_negbinomial)
               d = dzanbinom_v(y(row), parent, sizev, pzero)
            case (za_family_geometric)
               d = dzageom_v(y(row), parent, pzero)
            case default
               d = dzabinom_v(y(row), trials, parent, pzero)
            end select
            if (d <= tiny(1.0_dp) .or. .not. finite_scalar(d)) then
               nll = huge(1.0_dp)/100.0_dp; return
            end if
            nll = nll - w(row)*log(d)
         end do
      end function objective
   end subroutine fit_za_common

   subroutine predict_zero_altered_count(self, xp, xz, parent_parameter, zero_probability, fitted_mean)
      class(zero_altered_count_result_t), intent(in) :: self
      real(dp), intent(in) :: xp(:, :), xz(:, :)
      real(dp), allocatable, intent(out) :: parent_parameter(:), zero_probability(:), fitted_mean(:)
      integer :: n, i
      n = size(xp, 1)
      if (size(xz, 1) /= n .or. size(xp, 2) /= size(self%parent_coefficients) .or. &
          size(xz, 2) /= size(self%zero_coefficients)) then
         allocate(parent_parameter(0), zero_probability(0), fitted_mean(0)); return
      end if
      allocate(parent_parameter(n), zero_probability(n), fitted_mean(n))
      do i = 1, n
         call fitted_values_row(self, xp(i, :), xz(i, :), parent_parameter(i), zero_probability(i), fitted_mean(i))
      end do
   end subroutine predict_zero_altered_count

   subroutine fitted_values_row(self, xp, xz, parent, pzero, meanv)
      class(zero_altered_count_result_t), intent(in) :: self
      real(dp), intent(in) :: xp(:), xz(:)
      real(dp), intent(out) :: parent, pzero, meanv
      real(dp) :: base0, prob
      if (self%family == za_family_poisson .or. self%family == za_family_negbinomial) then
         parent = exp(clamp_eta(dot_product(xp, self%parent_coefficients)))
      else
         parent = logistic(dot_product(xp, self%parent_coefficients))
      end if
      pzero = logistic(dot_product(xz, self%zero_coefficients))
      select case (self%family)
      case (za_family_poisson)
         base0 = exp(-parent); meanv = (1.0_dp - pzero)*parent/max(1.0_dp - base0, tiny(1.0_dp))
      case (za_family_negbinomial)
         prob = self%size/(self%size + parent); base0 = prob**self%size
         meanv = (1.0_dp - pzero)*parent/max(1.0_dp - base0, tiny(1.0_dp))
      case (za_family_geometric)
         meanv = (1.0_dp - pzero)/max(parent, tiny(1.0_dp))
      case default
         base0 = (1.0_dp - parent)**self%trials
         meanv = (1.0_dp - pzero)*real(self%trials, dp)*parent/max(1.0_dp - base0, tiny(1.0_dp))
      end select
   end subroutine fitted_values_row

   elemental real(dp) function logistic(x) result(p)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then; p = 1.0_dp/(1.0_dp + exp(-x)); else; p = exp(x)/(1.0_dp + exp(x)); end if
   end function logistic

   elemental real(dp) function logit(p) result(x)
      real(dp), intent(in) :: p
      x = log(p/(1.0_dp - p))
   end function logit

   elemental real(dp) function clamp_eta(x) result(y)
      real(dp), intent(in) :: x
      y = min(30.0_dp, max(-30.0_dp, x))
   end function clamp_eta

   logical function is_intercept(x) result(ok)
      real(dp), intent(in) :: x(:, :)
      ok = size(x, 2) >= 1 .and. all(abs(x(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))
   end function is_intercept

   elemental logical function finite_scalar(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x) <= huge(x)
   end function finite_scalar
end module vgam_zero_altered_models
