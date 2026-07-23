! Part of the experimental modern Fortran translation of fGarch 4052.93.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original fGarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

module fgarch_distributions
   use fgarch_kinds, only : dp
   use fgarch_math, only : pi, normal_pdf, normal_cdf, normal_quantile, &
      student_t_pdf, student_t_cdf, student_t_quantile, regularized_gamma_p, &
      inverse_regularized_gamma_p, beta_fn
   use fgarch_rng, only : random_normal, random_gamma, random_student_t
   implicit none
   private

   integer, parameter, public :: dist_norm = 10
   integer, parameter, public :: dist_snorm = 11
   integer, parameter, public :: dist_std = 20
   integer, parameter, public :: dist_sstd = 21
   integer, parameter, public :: dist_ged = 30
   integer, parameter, public :: dist_sged = 31

   public :: dnorm_fg, pnorm_fg, qnorm_fg, rnorm_fg
   public :: dstd, pstd, qstd, rstd
   public :: dged, pged, qged, rged
   public :: dsnorm, psnorm, qsnorm, rsnorm
   public :: dsstd, psstd, qsstd, rsstd
   public :: dsged, psged, qsged, rsged
   public :: distribution_pdf, distribution_cdf, distribution_quantile
   public :: random_innovation, absolute_moment, distribution_name

contains

   pure elemental function dnorm_fg(x, mean, sd, log_density) result(value)
      real(dp), intent(in) :: x, mean, sd
      logical, intent(in), optional :: log_density
      real(dp) :: value
      logical :: want_log

      want_log = .false.
      if (present(log_density)) want_log = log_density
      if (sd <= 0.0_dp) then
         value = -huge(1.0_dp)
         if (.not. want_log) value = 0.0_dp
      else
         value = -0.5_dp*((x-mean)/sd)**2-log(sd)-0.5_dp*log(2.0_dp*pi)
         if (.not. want_log) value = exp(value)
      end if
   end function dnorm_fg

   pure elemental function pnorm_fg(q, mean, sd) result(value)
      real(dp), intent(in) :: q, mean, sd
      real(dp) :: value

      if (sd <= 0.0_dp) then
         value = 0.0_dp
      else
         value = normal_cdf((q-mean)/sd)
      end if
   end function pnorm_fg

   pure elemental function qnorm_fg(p, mean, sd) result(value)
      real(dp), intent(in) :: p, mean, sd
      real(dp) :: value

      value = mean + sd*normal_quantile(p)
   end function qnorm_fg

   function rnorm_fg(mean, sd) result(value)
      real(dp), intent(in) :: mean, sd
      real(dp) :: value

      value = mean + sd*random_normal()
   end function rnorm_fg

   pure elemental function dstd(x, mean, sd, nu, log_density) result(value)
      real(dp), intent(in) :: x, mean, sd, nu
      logical, intent(in), optional :: log_density
      real(dp) :: value, z, scale
      logical :: want_log

      want_log = .false.
      if (present(log_density)) want_log = log_density
      if (sd <= 0.0_dp .or. nu <= 2.0_dp) then
         value = -huge(1.0_dp)
         if (.not. want_log) value = 0.0_dp
         return
      end if
      scale = sqrt(nu/(nu-2.0_dp))
      z = (x-mean)/sd
      value = log(scale/sd) + log(student_t_pdf(z*scale,nu))
      if (.not. want_log) value = exp(value)
   end function dstd

   pure elemental function pstd(q, mean, sd, nu) result(value)
      real(dp), intent(in) :: q, mean, sd, nu
      real(dp) :: value, scale

      if (sd <= 0.0_dp .or. nu <= 2.0_dp) then
         value = 0.0_dp
      else
         scale = sqrt(nu/(nu-2.0_dp))
         value = student_t_cdf((q-mean)*scale/sd,nu)
      end if
   end function pstd

   pure elemental function qstd(p, mean, sd, nu) result(value)
      real(dp), intent(in) :: p, mean, sd, nu
      real(dp) :: value, scale

      scale = sqrt(nu/(nu-2.0_dp))
      value = mean + sd*student_t_quantile(p,nu)/scale
   end function qstd

   function rstd(mean, sd, nu) result(value)
      real(dp), intent(in) :: mean, sd, nu
      real(dp) :: value, scale

      scale = sqrt(nu/(nu-2.0_dp))
      value = mean + sd*random_student_t(nu)/scale
   end function rstd

   pure elemental function ged_lambda(nu) result(value)
      real(dp), intent(in) :: nu
      real(dp) :: value

      value = sqrt(2.0_dp**(-2.0_dp/nu)*gamma(1.0_dp/nu)/gamma(3.0_dp/nu))
   end function ged_lambda

   pure elemental function dged(x, mean, sd, nu, log_density) result(value)
      real(dp), intent(in) :: x, mean, sd, nu
      logical, intent(in), optional :: log_density
      real(dp) :: value, z, lambda, log_g
      logical :: want_log

      want_log = .false.
      if (present(log_density)) want_log = log_density
      if (sd <= 0.0_dp .or. nu <= 0.0_dp) then
         value = -huge(1.0_dp)
         if (.not. want_log) value = 0.0_dp
         return
      end if
      z = (x-mean)/sd
      lambda = ged_lambda(nu)
      log_g = log(nu)-log(lambda)-(1.0_dp+1.0_dp/nu)*log(2.0_dp)-log_gamma(1.0_dp/nu)
      value = log_g-0.5_dp*abs(z/lambda)**nu-log(sd)
      if (.not. want_log) value = exp(value)
   end function dged

   pure elemental function pged(q, mean, sd, nu) result(value)
      real(dp), intent(in) :: q, mean, sd, nu
      real(dp) :: value, z, lambda, pg

      if (sd <= 0.0_dp .or. nu <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      z = (q-mean)/sd
      if (abs(z) <= tiny(1.0_dp)) then
         value = 0.5_dp
      else
         lambda = ged_lambda(nu)
         pg = regularized_gamma_p(1.0_dp/nu,0.5_dp*abs(z/lambda)**nu)
         value = 0.5_dp + 0.5_dp*sign(1.0_dp,z)*pg
      end if
   end function pged

   pure elemental function qged(p, mean, sd, nu) result(value)
      real(dp), intent(in) :: p, mean, sd, nu
      real(dp) :: value, lambda, xg

      if (p <= 0.0_dp) then
         value = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         value = huge(1.0_dp)
      else if (abs(p-0.5_dp) <= epsilon(p)) then
         value = mean
      else
         lambda = ged_lambda(nu)
         xg = inverse_regularized_gamma_p(1.0_dp/nu,abs(2.0_dp*p-1.0_dp))
         value = mean + sd*sign(1.0_dp,p-0.5_dp)*lambda*(2.0_dp*xg)**(1.0_dp/nu)
      end if
   end function qged

   function rged(mean, sd, nu) result(value)
      real(dp), intent(in) :: mean, sd, nu
      real(dp) :: value, u, z

      call random_number(u)
      z = ged_lambda(nu)*(2.0_dp*random_gamma(1.0_dp/nu))**(1.0_dp/nu)
      if (u < 0.5_dp) z = -z
      value = mean + sd*z
   end function rged

   pure elemental function symmetric_abs_mean(kind, shape) result(m1)
      integer, intent(in) :: kind
      real(dp), intent(in) :: shape
      real(dp) :: m1

      select case (kind)
      case (dist_norm,dist_snorm)
         m1 = sqrt(2.0_dp/pi)
      case (dist_std,dist_sstd)
         m1 = 2.0_dp*sqrt(shape-2.0_dp)/(shape-1.0_dp)/beta_fn(0.5_dp,0.5_dp*shape)
      case (dist_ged,dist_sged)
         m1 = 2.0_dp**(1.0_dp/shape)*ged_lambda(shape)*gamma(2.0_dp/shape)/gamma(1.0_dp/shape)
      case default
         m1 = sqrt(2.0_dp/pi)
      end select
   end function symmetric_abs_mean

   pure elemental subroutine skew_standardization(kind, shape, xi, mu, sigma)
      integer, intent(in) :: kind
      real(dp), intent(in) :: shape, xi
      real(dp), intent(out) :: mu, sigma
      real(dp) :: m1

      m1 = symmetric_abs_mean(kind,shape)
      mu = m1*(xi-1.0_dp/xi)
      sigma = sqrt((1.0_dp-m1*m1)*(xi*xi+1.0_dp/(xi*xi))+2.0_dp*m1*m1-1.0_dp)
   end subroutine skew_standardization

   pure elemental function base_pdf(x, kind, shape) result(value)
      real(dp), intent(in) :: x, shape
      integer, intent(in) :: kind
      real(dp) :: value

      select case (kind)
      case (dist_norm,dist_snorm)
         value = normal_pdf(x)
      case (dist_std,dist_sstd)
         value = dstd(x,0.0_dp,1.0_dp,shape)
      case (dist_ged,dist_sged)
         value = dged(x,0.0_dp,1.0_dp,shape)
      case default
         value = normal_pdf(x)
      end select
   end function base_pdf

   pure elemental function base_cdf(x, kind, shape) result(value)
      real(dp), intent(in) :: x, shape
      integer, intent(in) :: kind
      real(dp) :: value

      select case (kind)
      case (dist_norm,dist_snorm)
         value = normal_cdf(x)
      case (dist_std,dist_sstd)
         value = pstd(x,0.0_dp,1.0_dp,shape)
      case (dist_ged,dist_sged)
         value = pged(x,0.0_dp,1.0_dp,shape)
      case default
         value = normal_cdf(x)
      end select
   end function base_cdf

   pure elemental function base_quantile(p, kind, shape) result(value)
      real(dp), intent(in) :: p, shape
      integer, intent(in) :: kind
      real(dp) :: value

      select case (kind)
      case (dist_norm,dist_snorm)
         value = normal_quantile(p)
      case (dist_std,dist_sstd)
         value = qstd(p,0.0_dp,1.0_dp,shape)
      case (dist_ged,dist_sged)
         value = qged(p,0.0_dp,1.0_dp,shape)
      case default
         value = normal_quantile(p)
      end select
   end function base_quantile

   pure elemental function skew_pdf_standard(x, kind, shape, xi) result(value)
      real(dp), intent(in) :: x, shape, xi
      integer, intent(in) :: kind
      real(dp) :: value, mu, sigma, z, xi_z, g

      if (xi <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      call skew_standardization(kind,shape,xi,mu,sigma)
      z = x*sigma+mu
      if (z > 0.0_dp) then
         xi_z = xi
      else if (z < 0.0_dp) then
         xi_z = 1.0_dp/xi
      else
         xi_z = 1.0_dp
      end if
      g = 2.0_dp/(xi+1.0_dp/xi)
      value = g*base_pdf(z/xi_z,kind,shape)*sigma
   end function skew_pdf_standard

   pure elemental function skew_cdf_standard(x, kind, shape, xi) result(value)
      real(dp), intent(in) :: x, shape, xi
      integer, intent(in) :: kind
      real(dp) :: value, mu, sigma, z, xi_z, g, sig

      call skew_standardization(kind,shape,xi,mu,sigma)
      z = x*sigma+mu
      if (z >= 0.0_dp) then
         sig = 1.0_dp
         xi_z = xi
         value = 1.0_dp
      else
         sig = -1.0_dp
         xi_z = 1.0_dp/xi
         value = 0.0_dp
      end if
      g = 2.0_dp/(xi+1.0_dp/xi)
      value = value-sig*g*xi_z*base_cdf(-abs(z)/xi_z,kind,shape)
      value = max(0.0_dp,min(1.0_dp,value))
   end function skew_cdf_standard

   pure elemental function skew_quantile_standard(p, kind, shape, xi) result(value)
      real(dp), intent(in) :: p, shape, xi
      integer, intent(in) :: kind
      real(dp) :: value, mu, sigma, g, p0, sig, xi_z, pp, h

      if (p <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      else if (p >= 1.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      call skew_standardization(kind,shape,xi,mu,sigma)
      p0 = 1.0_dp/(1.0_dp+xi*xi)
      if (abs(p-p0) <= 8.0_dp*epsilon(p)) then
         value = -mu/sigma
         return
      end if
      sig = sign(1.0_dp,p-p0)
      if (sig > 0.0_dp) then
         xi_z = xi
         h = 1.0_dp
      else
         xi_z = 1.0_dp/xi
         h = 0.0_dp
      end if
      g = 2.0_dp/(xi+1.0_dp/xi)
      pp = (h-sig*p)/(g*xi_z)
      pp = max(1.0e-15_dp,min(1.0_dp-1.0e-15_dp,pp))
      value = (-sig*xi_z*base_quantile(pp,kind,shape)-mu)/sigma
   end function skew_quantile_standard

   pure elemental function dsnorm(x, mean, sd, xi, log_density) result(value)
      real(dp), intent(in) :: x, mean, sd, xi
      logical, intent(in), optional :: log_density
      real(dp) :: value
      logical :: want_log

      want_log = .false.
      if (present(log_density)) want_log = log_density
      value = skew_pdf_standard((x-mean)/sd,dist_snorm,2.0_dp,xi)/sd
      if (want_log) value = log(max(value,tiny(1.0_dp)))
   end function dsnorm

   pure elemental function psnorm(q, mean, sd, xi) result(value)
      real(dp), intent(in) :: q, mean, sd, xi
      real(dp) :: value
      value = skew_cdf_standard((q-mean)/sd,dist_snorm,2.0_dp,xi)
   end function psnorm

   pure elemental function qsnorm(p, mean, sd, xi) result(value)
      real(dp), intent(in) :: p, mean, sd, xi
      real(dp) :: value
      value = mean+sd*skew_quantile_standard(p,dist_snorm,2.0_dp,xi)
   end function qsnorm

   function rsnorm(mean, sd, xi) result(value)
      real(dp), intent(in) :: mean, sd, xi
      real(dp) :: value
      value = mean+sd*random_skew_standard(dist_snorm,2.0_dp,xi)
   end function rsnorm

   pure elemental function dsstd(x, mean, sd, nu, xi, log_density) result(value)
      real(dp), intent(in) :: x, mean, sd, nu, xi
      logical, intent(in), optional :: log_density
      real(dp) :: value
      logical :: want_log
      want_log = .false.
      if (present(log_density)) want_log = log_density
      value = skew_pdf_standard((x-mean)/sd,dist_sstd,nu,xi)/sd
      if (want_log) value = log(max(value,tiny(1.0_dp)))
   end function dsstd

   pure elemental function psstd(q, mean, sd, nu, xi) result(value)
      real(dp), intent(in) :: q, mean, sd, nu, xi
      real(dp) :: value
      value = skew_cdf_standard((q-mean)/sd,dist_sstd,nu,xi)
   end function psstd

   pure elemental function qsstd(p, mean, sd, nu, xi) result(value)
      real(dp), intent(in) :: p, mean, sd, nu, xi
      real(dp) :: value
      value = mean+sd*skew_quantile_standard(p,dist_sstd,nu,xi)
   end function qsstd

   function rsstd(mean, sd, nu, xi) result(value)
      real(dp), intent(in) :: mean, sd, nu, xi
      real(dp) :: value
      value = mean+sd*random_skew_standard(dist_sstd,nu,xi)
   end function rsstd

   pure elemental function dsged(x, mean, sd, nu, xi, log_density) result(value)
      real(dp), intent(in) :: x, mean, sd, nu, xi
      logical, intent(in), optional :: log_density
      real(dp) :: value
      logical :: want_log
      want_log = .false.
      if (present(log_density)) want_log = log_density
      value = skew_pdf_standard((x-mean)/sd,dist_sged,nu,xi)/sd
      if (want_log) value = log(max(value,tiny(1.0_dp)))
   end function dsged

   pure elemental function psged(q, mean, sd, nu, xi) result(value)
      real(dp), intent(in) :: q, mean, sd, nu, xi
      real(dp) :: value
      value = skew_cdf_standard((q-mean)/sd,dist_sged,nu,xi)
   end function psged

   pure elemental function qsged(p, mean, sd, nu, xi) result(value)
      real(dp), intent(in) :: p, mean, sd, nu, xi
      real(dp) :: value
      value = mean+sd*skew_quantile_standard(p,dist_sged,nu,xi)
   end function qsged

   function rsged(mean, sd, nu, xi) result(value)
      real(dp), intent(in) :: mean, sd, nu, xi
      real(dp) :: value
      value = mean+sd*random_skew_standard(dist_sged,nu,xi)
   end function rsged

   function random_base(kind, shape) result(value)
      integer, intent(in) :: kind
      real(dp), intent(in) :: shape
      real(dp) :: value

      select case (kind)
      case (dist_norm,dist_snorm)
         value = random_normal()
      case (dist_std,dist_sstd)
         value = rstd(0.0_dp,1.0_dp,shape)
      case (dist_ged,dist_sged)
         value = rged(0.0_dp,1.0_dp,shape)
      case default
         value = random_normal()
      end select
   end function random_base

   function random_skew_standard(kind, shape, xi) result(value)
      integer, intent(in) :: kind
      real(dp), intent(in) :: shape, xi
      real(dp) :: value, u, weight, sig, xi_z, mu, sigma

      weight = xi/(xi+1.0_dp/xi)
      call random_number(u)
      u = u-weight
      if (u >= 0.0_dp) then
         sig = 1.0_dp
         xi_z = xi
      else
         sig = -1.0_dp
         xi_z = 1.0_dp/xi
      end if
      value = -abs(random_base(kind,shape))/xi_z*sig
      call skew_standardization(kind,shape,xi,mu,sigma)
      value = (value-mu)/sigma
   end function random_skew_standard

   pure elemental function distribution_pdf(x, kind, shape, skew) result(value)
      real(dp), intent(in) :: x, shape, skew
      integer, intent(in) :: kind
      real(dp) :: value

      select case (kind)
      case (dist_norm)
         value = normal_pdf(x)
      case (dist_snorm)
         value = dsnorm(x,0.0_dp,1.0_dp,skew)
      case (dist_std)
         value = dstd(x,0.0_dp,1.0_dp,shape)
      case (dist_sstd)
         value = dsstd(x,0.0_dp,1.0_dp,shape,skew)
      case (dist_ged)
         value = dged(x,0.0_dp,1.0_dp,shape)
      case (dist_sged)
         value = dsged(x,0.0_dp,1.0_dp,shape,skew)
      case default
         value = normal_pdf(x)
      end select
   end function distribution_pdf

   pure elemental function distribution_cdf(x, kind, shape, skew) result(value)
      real(dp), intent(in) :: x, shape, skew
      integer, intent(in) :: kind
      real(dp) :: value

      select case (kind)
      case (dist_norm)
         value = normal_cdf(x)
      case (dist_snorm)
         value = psnorm(x,0.0_dp,1.0_dp,skew)
      case (dist_std)
         value = pstd(x,0.0_dp,1.0_dp,shape)
      case (dist_sstd)
         value = psstd(x,0.0_dp,1.0_dp,shape,skew)
      case (dist_ged)
         value = pged(x,0.0_dp,1.0_dp,shape)
      case (dist_sged)
         value = psged(x,0.0_dp,1.0_dp,shape,skew)
      case default
         value = normal_cdf(x)
      end select
   end function distribution_cdf

   pure elemental function distribution_quantile(p, kind, shape, skew) result(value)
      real(dp), intent(in) :: p, shape, skew
      integer, intent(in) :: kind
      real(dp) :: value

      select case (kind)
      case (dist_norm)
         value = normal_quantile(p)
      case (dist_snorm)
         value = qsnorm(p,0.0_dp,1.0_dp,skew)
      case (dist_std)
         value = qstd(p,0.0_dp,1.0_dp,shape)
      case (dist_sstd)
         value = qsstd(p,0.0_dp,1.0_dp,shape,skew)
      case (dist_ged)
         value = qged(p,0.0_dp,1.0_dp,shape)
      case (dist_sged)
         value = qsged(p,0.0_dp,1.0_dp,shape,skew)
      case default
         value = normal_quantile(p)
      end select
   end function distribution_quantile

   function random_innovation(kind, shape, skew) result(value)
      integer, intent(in) :: kind
      real(dp), intent(in) :: shape, skew
      real(dp) :: value

      select case (kind)
      case (dist_norm)
         value = random_normal()
      case (dist_snorm)
         value = rsnorm(0.0_dp,1.0_dp,skew)
      case (dist_std)
         value = rstd(0.0_dp,1.0_dp,shape)
      case (dist_sstd)
         value = rsstd(0.0_dp,1.0_dp,shape,skew)
      case (dist_ged)
         value = rged(0.0_dp,1.0_dp,shape)
      case (dist_sged)
         value = rsged(0.0_dp,1.0_dp,shape,skew)
      case default
         value = random_normal()
      end select
   end function random_innovation

   pure elemental function absolute_moment(order, kind, shape) result(value)
      real(dp), intent(in) :: order, shape
      integer, intent(in) :: kind
      real(dp) :: value

      select case (kind)
      case (dist_norm,dist_snorm)
         value = sqrt(2.0_dp)**order*gamma(0.5_dp*(order+1.0_dp))/sqrt(pi)
      case (dist_ged,dist_sged)
         value = (2.0_dp**(1.0_dp/shape)*ged_lambda(shape))**order * &
                 gamma((order+1.0_dp)/shape)/gamma(1.0_dp/shape)
      case (dist_std,dist_sstd)
         if (shape <= order) then
            value = huge(1.0_dp)
         else
            value = beta_fn(0.5_dp+0.5_dp*order,0.5_dp*shape-0.5_dp*order) / &
                    beta_fn(0.5_dp,0.5_dp*shape)*(shape-2.0_dp)**(0.5_dp*order)
         end if
      case default
         value = 0.0_dp
      end select
   end function absolute_moment

   pure function distribution_name(kind) result(name)
      integer, intent(in) :: kind
      character(len=8) :: name

      select case (kind)
      case (dist_norm);  name = 'norm'
      case (dist_snorm); name = 'snorm'
      case (dist_std);   name = 'std'
      case (dist_sstd);  name = 'sstd'
      case (dist_ged);   name = 'ged'
      case (dist_sged);  name = 'sged'
      case default;      name = 'unknown'
      end select
   end function distribution_name

end module fgarch_distributions
