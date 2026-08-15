! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_copulas
   use vgam_kinds, only : dp
   use vgam_distributions, only : qnorm_v, pnorm_v, rnorm_v
   use vgam_special, only : log1p_v
   use vgam_optim, only : bfgs_minimize, numerical_hessian
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   integer, parameter, public :: copula_clayton = 1
   integer, parameter, public :: copula_frank = 2
   integer, parameter, public :: copula_fgm = 3
   integer, parameter, public :: copula_gaussian = 4
   integer, parameter, public :: copula_plackett = 5
   integer, parameter, public :: copula_amh = 6

   type, public :: copula_regression_result_t
      integer :: family = 0
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted_parameter(:)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: status = 0
      logical :: converged = .false.
   contains
      procedure :: predict_parameter => predict_copula_parameter
   end type copula_regression_result_t

   public :: clayton_copula_pdf, clayton_copula_cdf, random_clayton_copula
   public :: frank_copula_pdf, frank_copula_cdf, random_frank_copula
   public :: fgm_copula_pdf, fgm_copula_cdf, random_fgm_copula
   public :: gaussian_copula_pdf, gaussian_copula_cdf, random_gaussian_copula
   public :: plackett_copula_pdf, plackett_copula_cdf, random_plackett_copula
   public :: amh_copula_pdf, amh_copula_cdf, random_amh_copula
   public :: fit_copula_regression

contains

   elemental real(dp) function clayton_copula_cdf(u, v, apar) result(c)
      real(dp), intent(in) :: u, v, apar
      real(dp) :: a, s
      if (u <= 0.0_dp .or. v <= 0.0_dp) then
         c = 0.0_dp
      else if (u >= 1.0_dp .and. v >= 1.0_dp) then
         c = 1.0_dp
      else if (u >= 1.0_dp) then
         c = min(1.0_dp, v)
      else if (v >= 1.0_dp) then
         c = min(1.0_dp, u)
      else if (apar < 0.0_dp) then
         c = ieee_nan()
      else if (apar <= sqrt(epsilon(1.0_dp))) then
         c = u*v
      else
         a = apar
         s = u**(-a) + v**(-a) - 1.0_dp
         if (s <= 0.0_dp) then
            c = ieee_nan()
         else
            c = s**(-1.0_dp/a)
         end if
      end if
   end function clayton_copula_cdf

   elemental real(dp) function clayton_copula_pdf(u, v, apar) result(d)
      real(dp), intent(in) :: u, v, apar
      real(dp) :: a, s, ld
      if (u <= 0.0_dp .or. u >= 1.0_dp .or. v <= 0.0_dp .or. v >= 1.0_dp) then
         d = 0.0_dp
      else if (apar < 0.0_dp) then
         d = ieee_nan()
      else if (apar <= sqrt(epsilon(1.0_dp))) then
         d = 1.0_dp
      else
         a = apar
         s = u**(-a) + v**(-a) - 1.0_dp
         if (s <= 0.0_dp) then
            d = ieee_nan()
         else
            ld = log1p_v(a) - (1.0_dp + a)*(log(u) + log(v)) &
                 - (2.0_dp + 1.0_dp/a)*log(s)
            d = exp(min(700.0_dp, ld))
         end if
      end if
   end function clayton_copula_pdf

   subroutine random_clayton_copula(apar, u, v)
      real(dp), intent(in) :: apar
      real(dp), intent(out) :: u, v
      real(dp) :: w
      call random_number(u)
      call random_number(w)
      u = interior_uniform(u)
      w = interior_uniform(w)
      if (apar <= sqrt(epsilon(1.0_dp))) then
         v = w
      else
         v = (u**(-apar)*(w**(-apar/(1.0_dp + apar)) - 1.0_dp) + 1.0_dp)**(-1.0_dp/apar)
      end if
   end subroutine random_clayton_copula

   elemental real(dp) function frank_copula_cdf(u, v, apar) result(c)
      real(dp), intent(in) :: u, v, apar
      real(dp) :: la, term
      if (u <= 0.0_dp .or. v <= 0.0_dp) then
         c = 0.0_dp
      else if (u >= 1.0_dp .and. v >= 1.0_dp) then
         c = 1.0_dp
      else if (u >= 1.0_dp) then
         c = min(1.0_dp, v)
      else if (v >= 1.0_dp) then
         c = min(1.0_dp, u)
      else if (apar <= 0.0_dp) then
         c = ieee_nan()
      else if (abs(apar - 1.0_dp) <= 20.0_dp*epsilon(1.0_dp)) then
         c = u*v
      else
         la = log(apar)
         term = 1.0_dp + (exp(la*u) - 1.0_dp)*(exp(la*v) - 1.0_dp)/(apar - 1.0_dp)
         if (term <= 0.0_dp) then
            c = ieee_nan()
         else
            c = log(term)/la
         end if
      end if
   end function frank_copula_cdf

   elemental real(dp) function frank_copula_pdf(u, v, apar) result(d)
      real(dp), intent(in) :: u, v, apar
      real(dp) :: la, temp, num
      if (u <= 0.0_dp .or. u >= 1.0_dp .or. v <= 0.0_dp .or. v >= 1.0_dp) then
         d = 0.0_dp
      else if (apar <= 0.0_dp) then
         d = ieee_nan()
      else if (abs(apar - 1.0_dp) <= 20.0_dp*epsilon(1.0_dp)) then
         d = 1.0_dp
      else
         la = log(apar)
         temp = (apar - 1.0_dp) + (exp(la*u) - 1.0_dp)*(exp(la*v) - 1.0_dp)
         num = (apar - 1.0_dp)*la*exp(la*(u + v))
         d = num/(temp*temp)
      end if
   end function frank_copula_pdf

   subroutine random_frank_copula(apar, u, v)
      real(dp), intent(in) :: apar
      real(dp), intent(out) :: u, v
      real(dp) :: w, t, la
      call random_number(u)
      call random_number(w)
      u = interior_uniform(u)
      w = interior_uniform(w)
      if (apar <= 0.0_dp) then
         v = ieee_nan()
      else if (abs(apar - 1.0_dp) <= 20.0_dp*epsilon(1.0_dp)) then
         v = w
      else
         la = log(apar)
         t = exp(la*u) + (apar - exp(la*u))*w
         v = log(t/(t + (1.0_dp - apar)*w))/la
      end if
   end subroutine random_frank_copula

   elemental real(dp) function fgm_copula_cdf(u, v, apar) result(c)
      real(dp), intent(in) :: u, v, apar
      if (abs(apar) > 1.0_dp) then
         c = ieee_nan()
      else if (u <= 0.0_dp .or. v <= 0.0_dp) then
         c = 0.0_dp
      else if (u >= 1.0_dp .and. v >= 1.0_dp) then
         c = 1.0_dp
      else if (u >= 1.0_dp) then
         c = min(1.0_dp, v)
      else if (v >= 1.0_dp) then
         c = min(1.0_dp, u)
      else
         c = u*v*(1.0_dp + apar*(1.0_dp - u)*(1.0_dp - v))
      end if
   end function fgm_copula_cdf

   elemental real(dp) function fgm_copula_pdf(u, v, apar) result(d)
      real(dp), intent(in) :: u, v, apar
      if (abs(apar) > 1.0_dp) then
         d = ieee_nan()
      else if (u <= 0.0_dp .or. u >= 1.0_dp .or. v <= 0.0_dp .or. v >= 1.0_dp) then
         d = 0.0_dp
      else
         d = 1.0_dp + apar*(1.0_dp - 2.0_dp*u)*(1.0_dp - 2.0_dp*v)
      end if
   end function fgm_copula_pdf

   subroutine random_fgm_copula(apar, u, v)
      real(dp), intent(in) :: apar
      real(dp), intent(out) :: u, v
      real(dp) :: w, temp, aa, bb, disc
      call random_number(u)
      call random_number(w)
      u = interior_uniform(u)
      w = interior_uniform(w)
      if (abs(apar) <= sqrt(epsilon(1.0_dp))) then
         v = w
      else
         temp = 2.0_dp*u - 1.0_dp
         aa = apar*temp - 1.0_dp
         disc = 1.0_dp - 2.0_dp*apar*temp + (apar*temp)**2 + 4.0_dp*apar*w*temp
         bb = sqrt(max(0.0_dp, disc))
         v = 2.0_dp*w/(bb - aa)
      end if
   end subroutine random_fgm_copula

   elemental real(dp) function gaussian_copula_pdf(u, v, rho) result(d)
      real(dp), intent(in) :: u, v, rho
      real(dp) :: z1, z2, den, ld
      if (u <= 0.0_dp .or. u >= 1.0_dp .or. v <= 0.0_dp .or. v >= 1.0_dp) then
         d = 0.0_dp
      else if (abs(rho) >= 1.0_dp) then
         d = ieee_nan()
      else
         z1 = qnorm_v(u, 0.0_dp, 1.0_dp)
         z2 = qnorm_v(v, 0.0_dp, 1.0_dp)
         den = 1.0_dp - rho*rho
         ld = (2.0_dp*rho*z1*z2 - rho*rho*(z1*z1 + z2*z2))/(2.0_dp*den) &
              - 0.5_dp*log(den)
         d = exp(min(700.0_dp, ld))
      end if
   end function gaussian_copula_pdf

   real(dp) function gaussian_copula_cdf(u, v, rho) result(c)
      real(dp), intent(in) :: u, v, rho
      real(dp) :: a, b, lo, h, x, s, den
      integer, parameter :: nstep = 800
      integer :: i
      if (u <= 0.0_dp .or. v <= 0.0_dp) then
         c = 0.0_dp
         return
      else if (u >= 1.0_dp .and. v >= 1.0_dp) then
         c = 1.0_dp
         return
      else if (u >= 1.0_dp) then
         c = min(1.0_dp, v)
         return
      else if (v >= 1.0_dp) then
         c = min(1.0_dp, u)
         return
      else if (abs(rho) >= 1.0_dp) then
         c = ieee_nan()
         return
      end if
      if (abs(rho) <= 10.0_dp*epsilon(1.0_dp)) then
         c = u*v
         return
      end if
      a = qnorm_v(u, 0.0_dp, 1.0_dp)
      b = qnorm_v(v, 0.0_dp, 1.0_dp)
      lo = -9.0_dp
      if (a <= lo) then
         c = 0.0_dp
         return
      end if
      den = sqrt(1.0_dp - rho*rho)
      h = (a - lo)/real(nstep, dp)
      s = integrand(lo, b, rho, den) + integrand(a, b, rho, den)
      do i = 1, nstep - 1
         x = lo + real(i, dp)*h
         if (mod(i, 2) == 0) then
            s = s + 2.0_dp*integrand(x, b, rho, den)
         else
            s = s + 4.0_dp*integrand(x, b, rho, den)
         end if
      end do
      c = min(1.0_dp, max(0.0_dp, h*s/3.0_dp))
   end function gaussian_copula_cdf

   subroutine random_gaussian_copula(rho, u, v)
      real(dp), intent(in) :: rho
      real(dp), intent(out) :: u, v
      real(dp) :: z1, z2
      if (abs(rho) >= 1.0_dp) then
         u = ieee_nan()
         v = ieee_nan()
         return
      end if
      z1 = rnorm_v(0.0_dp, 1.0_dp)
      z2 = rho*z1 + sqrt(1.0_dp - rho*rho)*rnorm_v(0.0_dp, 1.0_dp)
      u = pnorm_v(z1, 0.0_dp, 1.0_dp)
      v = pnorm_v(z2, 0.0_dp, 1.0_dp)
      u = interior_uniform(u)
      v = interior_uniform(v)
   end subroutine random_gaussian_copula


   elemental real(dp) function plackett_copula_cdf(u, v, oratio) result(c)
      real(dp), intent(in) :: u, v, oratio
      real(dp) :: temp1, disc
      if (oratio <= 0.0_dp) then
         c = ieee_nan()
      else if (u <= 0.0_dp .or. v <= 0.0_dp) then
         c = 0.0_dp
      else if (u >= 1.0_dp .and. v >= 1.0_dp) then
         c = 1.0_dp
      else if (u >= 1.0_dp) then
         c = min(1.0_dp, v)
      else if (v >= 1.0_dp) then
         c = min(1.0_dp, u)
      else if (abs(oratio - 1.0_dp) < 1.0e-6_dp) then
         c = u*v
      else
         temp1 = 1.0_dp + (oratio - 1.0_dp)*(u + v)
         disc = temp1*temp1 - 4.0_dp*oratio*(oratio - 1.0_dp)*u*v
         c = 0.5_dp*(temp1 - sqrt(max(0.0_dp, disc)))/(oratio - 1.0_dp)
      end if
   end function plackett_copula_cdf

   elemental real(dp) function plackett_copula_pdf(u, v, oratio) result(d)
      real(dp), intent(in) :: u, v, oratio
      real(dp) :: num, den, ld
      if (oratio <= 0.0_dp) then
         d = ieee_nan()
      else if (u <= 0.0_dp .or. u >= 1.0_dp .or. v <= 0.0_dp .or. v >= 1.0_dp) then
         d = 0.0_dp
      else if (abs(oratio - 1.0_dp) < 1.0e-8_dp) then
         d = 1.0_dp
      else
         num = 1.0_dp + (oratio - 1.0_dp)*(u + v - 2.0_dp*u*v)
         den = (1.0_dp + (u + v)*(oratio - 1.0_dp))**2 &
               - 4.0_dp*oratio*(oratio - 1.0_dp)*u*v
         if (den <= 0.0_dp .or. num <= 0.0_dp) then
            d = ieee_nan()
         else
            ld = log(oratio) + log(num) - 1.5_dp*log(den)
            d = exp(min(700.0_dp, ld))
         end if
      end if
   end function plackett_copula_pdf

   subroutine random_plackett_copula(oratio, u, v)
      real(dp), intent(in) :: oratio
      real(dp), intent(out) :: u, v
      real(dp) :: w, z, y2, rad
      call random_number(u)
      call random_number(w)
      u = interior_uniform(u)
      w = interior_uniform(w)
      if (oratio <= 0.0_dp) then
         v = ieee_nan()
         return
      end if
      if (abs(oratio - 1.0_dp) < 1.0e-8_dp) then
         v = w
         return
      end if
      z = w*(1.0_dp - w)
      rad = oratio*(oratio + 4.0_dp*z*u*(1.0_dp - u)*(1.0_dp - oratio)**2)
      y2 = (2.0_dp*z*(u*oratio*oratio + 1.0_dp - u) + oratio*(1.0_dp - 2.0_dp*z) &
            - (1.0_dp - 2.0_dp*w)*sqrt(max(0.0_dp, rad))) &
            /(oratio + z*(1.0_dp - oratio)**2)
      v = 0.5_dp*y2
      v = interior_uniform(v)
   end subroutine random_plackett_copula


   elemental real(dp) function amh_copula_cdf(u, v, apar) result(c)
      real(dp), intent(in) :: u, v, apar
      real(dp) :: den
      if (abs(apar) > 1.0_dp) then
         c = ieee_nan()
      else if (u <= 0.0_dp .or. v <= 0.0_dp) then
         c = 0.0_dp
      else if (u >= 1.0_dp .and. v >= 1.0_dp) then
         c = 1.0_dp
      else if (u >= 1.0_dp) then
         c = min(1.0_dp, v)
      else if (v >= 1.0_dp) then
         c = min(1.0_dp, u)
      else
         den = 1.0_dp - apar*(1.0_dp - u)*(1.0_dp - v)
         c = u*v/den
      end if
   end function amh_copula_cdf

   elemental real(dp) function amh_copula_pdf(u, v, apar) result(d)
      real(dp), intent(in) :: u, v, apar
      real(dp) :: temp
      if (abs(apar) > 1.0_dp) then
         d = ieee_nan()
      else if (u <= 0.0_dp .or. u >= 1.0_dp .or. v <= 0.0_dp .or. v >= 1.0_dp) then
         d = 0.0_dp
      else
         temp = 1.0_dp - apar*(1.0_dp - u)*(1.0_dp - v)
         d = (1.0_dp - apar + 2.0_dp*apar*u*v/temp)/(temp*temp)
      end if
   end function amh_copula_pdf

   subroutine random_amh_copula(apar, u, v)
      real(dp), intent(in) :: apar
      real(dp), intent(out) :: u, v
      real(dp) :: v1, v2, b, aa, bb, rad
      call random_number(v1)
      call random_number(v2)
      u = interior_uniform(v1)
      if (abs(apar) > 1.0_dp) then
         v = ieee_nan()
         return
      else if (abs(apar) <= 20.0_dp*epsilon(1.0_dp)) then
         v = interior_uniform(v2)
         return
      end if
      b = 1.0_dp - v1
      aa = -apar*(2.0_dp*b*v2 + 1.0_dp) + 2.0_dp*apar*apar*b*b*v2 + 1.0_dp
      bb = apar*apar*(4.0_dp*b*b*v2 - 4.0_dp*b*v2 + 1.0_dp) &
           + apar*(4.0_dp*v2 - 4.0_dp*b*v2 - 2.0_dp) + 1.0_dp
      rad = max(0.0_dp, bb)
      v = 2.0_dp*v2*(apar*b - 1.0_dp)**2/(aa + sqrt(rad))
      v = interior_uniform(v)
   end subroutine random_amh_copula

   subroutine fit_copula_regression(u, v, x, family, result, weights, max_iter, tol)
      real(dp), intent(in) :: u(:), v(:), x(:, :)
      integer, intent(in) :: family
      type(copula_regression_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), tol
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: w(:), par(:), hess(:, :), cov(:, :)
      real(dp) :: fval, tolerance
      integer :: n, p, stat, stat2, niter, i

      n = size(u)
      p = size(x, 2)
      if (n <= 1 .or. size(v) /= n .or. size(x, 1) /= n .or. p <= 0 .or. &
          any(u <= 0.0_dp) .or. any(u >= 1.0_dp) .or. any(v <= 0.0_dp) .or. any(v >= 1.0_dp) .or. &
          family < copula_clayton .or. family > copula_amh) then
         result%status = 1
         return
      end if
      allocate(w(n))
      w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            result%status = 2
            return
         end if
         w = weights
      end if
      allocate(par(p))
      par = 0.0_dp
      if (all(abs(x(:, 1) - 1.0_dp) < 100.0_dp*epsilon(1.0_dp))) then
         select case (family)
         case (copula_clayton)
            par(1) = log(0.25_dp)
         case (copula_fgm)
            par(1) = atanh(clamp_dependence(3.0_dp*sample_correlation(u, v), 0.85_dp))
         case (copula_gaussian)
            par(1) = atanh(clamp_dependence(normal_score_correlation(u, v), 0.95_dp))
         case (copula_amh)
            par(1) = atanh(clamp_dependence(3.0_dp*sample_correlation(u, v), 0.85_dp))
         case default
            par(1) = 0.0_dp
         end select
      end if
      niter = 300
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp
      if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, par, fval, stat, max_iter=niter, tol=tolerance)
      result%family = family
      result%coefficients = par
      result%status = stat
      result%converged = stat == 0
      result%loglik = -fval
      result%aic = 2.0_dp*fval + 2.0_dp*real(p, dp)
      allocate(result%fitted_parameter(n))
      do i = 1, n
         result%fitted_parameter(i) = eta_to_parameter(dot_product(x(i, :), par), family)
      end do
      allocate(hess(p, p))
      call numerical_hessian(objective, par, hess)
      call invert_matrix(hess, cov, stat2)
      if (stat2 == 0) then
         result%covariance = cov
      else
         allocate(result%covariance(0, 0))
      end if

   contains

      real(dp) function objective(beta) result(nll)
         real(dp), intent(in) :: beta(:)
         real(dp) :: apar, dens
         integer :: j
         nll = 0.0_dp
         do j = 1, n
            apar = eta_to_parameter(dot_product(x(j, :), beta), family)
            dens = copula_density(u(j), v(j), apar, family)
            if (dens <= tiny(1.0_dp) .or. .not. finite_scalar(dens)) then
               nll = huge(1.0_dp)/100.0_dp
               return
            end if
            nll = nll - w(j)*log(dens)
         end do
      end function objective

   end subroutine fit_copula_regression

   subroutine predict_copula_parameter(self, x, parameter)
      class(copula_regression_result_t), intent(in) :: self
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: parameter(:)
      integer :: i
      if (.not. allocated(self%coefficients) .or. size(x, 2) /= size(self%coefficients)) then
         allocate(parameter(0))
         return
      end if
      allocate(parameter(size(x, 1)))
      do i = 1, size(x, 1)
         parameter(i) = eta_to_parameter(dot_product(x(i, :), self%coefficients), self%family)
      end do
   end subroutine predict_copula_parameter

   elemental real(dp) function copula_density(u, v, apar, family) result(d)
      real(dp), intent(in) :: u, v, apar
      integer, intent(in) :: family
      select case (family)
      case (copula_clayton)
         d = clayton_copula_pdf(u, v, apar)
      case (copula_frank)
         d = frank_copula_pdf(u, v, apar)
      case (copula_fgm)
         d = fgm_copula_pdf(u, v, apar)
      case (copula_gaussian)
         d = gaussian_copula_pdf(u, v, apar)
      case (copula_plackett)
         d = plackett_copula_pdf(u, v, apar)
      case (copula_amh)
         d = amh_copula_pdf(u, v, apar)
      case default
         d = ieee_nan()
      end select
   end function copula_density

   elemental real(dp) function eta_to_parameter(eta, family) result(ap)
      real(dp), intent(in) :: eta
      integer, intent(in) :: family
      select case (family)
      case (copula_clayton, copula_frank, copula_plackett)
         ap = exp(min(30.0_dp, max(-30.0_dp, eta)))
      case (copula_fgm, copula_gaussian, copula_amh)
         ap = tanh(eta)
      case default
         ap = ieee_nan()
      end select
   end function eta_to_parameter

   pure real(dp) function integrand(x, b, rho, den) result(v)
      real(dp), intent(in) :: x, b, rho, den
      real(dp), parameter :: invsqrt2pi = 0.398942280401432677939946059934_dp
      v = invsqrt2pi*exp(-0.5_dp*x*x)*pnorm_v((b - rho*x)/den, 0.0_dp, 1.0_dp)
   end function integrand

   real(dp) function sample_correlation(x, y) result(r)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: mx, my, sx, sy
      mx = sum(x)/real(size(x), dp)
      my = sum(y)/real(size(y), dp)
      sx = sum((x - mx)**2)
      sy = sum((y - my)**2)
      if (sx <= tiny(1.0_dp) .or. sy <= tiny(1.0_dp)) then
         r = 0.0_dp
      else
         r = sum((x - mx)*(y - my))/sqrt(sx*sy)
      end if
   end function sample_correlation

   real(dp) function normal_score_correlation(u, v) result(r)
      real(dp), intent(in) :: u(:), v(:)
      real(dp), allocatable :: z1(:), z2(:)
      integer :: i
      allocate(z1(size(u)), z2(size(v)))
      do i = 1, size(u)
         z1(i) = qnorm_v(interior_uniform(u(i)), 0.0_dp, 1.0_dp)
         z2(i) = qnorm_v(interior_uniform(v(i)), 0.0_dp, 1.0_dp)
      end do
      r = sample_correlation(z1, z2)
   end function normal_score_correlation

   elemental real(dp) function clamp_dependence(x, lim) result(y)
      real(dp), intent(in) :: x, lim
      y = min(lim, max(-lim, x))
   end function clamp_dependence

   elemental real(dp) function interior_uniform(x) result(y)
      real(dp), intent(in) :: x
      y = min(1.0_dp - 1.0e-14_dp, max(1.0e-14_dp, x))
   end function interior_uniform

   elemental logical function finite_scalar(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x) <= huge(x)
   end function finite_scalar

   elemental real(dp) function ieee_nan() result(x)
      real(dp) :: z
      z = 0.0_dp
      x = 0.0_dp/z
   end function ieee_nan

end module vgam_copulas
