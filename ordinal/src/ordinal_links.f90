! Based on ordinal/src/links.c and ordinal/R/{gdist,gumbel,lgamma,AO}.R
! Copyright (C) 2011-2026 R. H. B. Christensen
! Modern Fortran translation, 2026. Distributed under GPL-2.0-or-later.
module ordinal_links
   use ordinal_kinds, only : dp
   use ordinal_special, only : normal_cdf, normal_pdf, regularized_gamma_p
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   implicit none
   private
   integer, parameter, public :: link_logit = 1
   integer, parameter, public :: link_probit = 2
   integer, parameter, public :: link_cloglog = 3
   integer, parameter, public :: link_loglog = 4
   integer, parameter, public :: link_cauchit = 5
   integer, parameter, public :: link_aranda_ordaz = 6
   integer, parameter, public :: link_log_gamma = 7
   real(dp), parameter :: pi = acos(-1.0_dp)
   public :: link_cdf, link_pdf, link_pdf_gradient
   public :: pgumbel, dgumbel, qgumbel, glogis, gnorm, gcauchy, ggumbel
   public :: paranda_ordaz, daranda_ordaz, garanda_ordaz
   public :: plgamma_link, dlgamma_link, glgamma_link
contains
   pure elemental real(dp) function pgumbel(q, location, scale, lower_tail, maximum) result(p)
      real(dp), intent(in) :: q !! Quantile at which the Gumbel CDF is evaluated.
      real(dp), intent(in) :: location !! Gumbel location parameter.
      real(dp), intent(in) :: scale !! Positive Gumbel scale parameter.
      logical, intent(in) :: lower_tail !! Return the lower tail when true and the upper tail otherwise.
      logical, intent(in) :: maximum !! Use the maximum-type Gumbel when true and minimum-type when false.
      real(dp) :: z
      if (ieee_is_nan(q) .or. scale <= 0.0_dp) then
         p = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (maximum) then
         z = (q - location)/scale
         p = exp(-exp(-z))
      else
         z = (q - location)/scale
         p = 1.0_dp - exp(-exp(z))
      end if
      if (.not. lower_tail) p = 1.0_dp - p
   end function pgumbel

   pure elemental real(dp) function dgumbel(x, location, scale, give_log, maximum) result(y)
      real(dp), intent(in) :: x !! Point at which the Gumbel density is evaluated.
      real(dp), intent(in) :: location !! Gumbel location parameter.
      real(dp), intent(in) :: scale !! Positive Gumbel scale parameter.
      logical, intent(in) :: give_log !! Return the log-density when true.
      logical, intent(in) :: maximum !! Use the maximum-type Gumbel when true and minimum-type when false.
      real(dp) :: z, ly
      if (ieee_is_nan(x) .or. scale <= 0.0_dp) then
         y = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (maximum) then
         z = (x - location)/scale
         ly = -exp(-z) - z - log(scale)
      else
         z = -(x - location)/scale
         ly = -exp(-z) - z - log(scale)
      end if
      if (give_log) then
         y = ly
      else
         y = exp(ly)
      end if
   end function dgumbel

   pure elemental real(dp) function qgumbel(p, location, scale, lower_tail, maximum) result(q)
      real(dp), intent(in) :: p !! Probability in the closed unit interval.
      real(dp), intent(in) :: location !! Gumbel location parameter.
      real(dp), intent(in) :: scale !! Positive Gumbel scale parameter.
      logical, intent(in) :: lower_tail !! Interpret p as a lower-tail probability when true.
      logical, intent(in) :: maximum !! Use the maximum-type Gumbel when true and minimum-type when false.
      real(dp) :: pp
      pp = merge(p, 1.0_dp - p, lower_tail)
      if (maximum) then
         q = location - scale*log(-log(pp))
      else
         q = location + scale*log(-log(1.0_dp - pp))
      end if
   end function qgumbel

   pure elemental real(dp) function glogis(x) result(g)
      real(dp), intent(in) :: x !! Logistic variate at which the density derivative is evaluated.
      real(dp) :: e
      e = exp(-abs(x))
      g = 2.0_dp*e*e/(1.0_dp + e)**3 - e/(1.0_dp + e)**2
      if (x <= 0.0_dp) g = -g
   end function glogis

   pure elemental real(dp) function gnorm(x) result(g)
      real(dp), intent(in) :: x !! Standard-normal variate at which the density derivative is evaluated.
      g = -x*normal_pdf(x)
   end function gnorm

   pure elemental real(dp) function gcauchy(x) result(g)
      real(dp), intent(in) :: x !! Standard-Cauchy variate at which the density derivative is evaluated.
      g = -2.0_dp*x/(pi*(1.0_dp + x*x)**2)
   end function gcauchy

   pure elemental real(dp) function ggumbel(x, maximum) result(g)
      real(dp), intent(in) :: x !! Standard Gumbel variate at which the density derivative is evaluated.
      logical, intent(in) :: maximum !! Use the maximum-type Gumbel when true and minimum-type when false.
      real(dp) :: z, ez
      z = merge(x, -x, maximum)
      ez = exp(-z)
      g = exp(-ez)*ez*(ez - 1.0_dp)
      if (.not. maximum) g = -g
   end function ggumbel

   pure elemental real(dp) function paranda_ordaz(eta, lambda, lower_tail) result(p)
      real(dp), intent(in) :: eta !! Linear predictor for the Aranda-Ordaz CDF.
      real(dp), intent(in) :: lambda !! Positive Aranda-Ordaz shape parameter.
      logical, intent(in) :: lower_tail !! Return the lower tail when true and upper tail otherwise.
      if (lambda <= 0.0_dp) then
         p = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      p = 1.0_dp - (lambda*exp(eta) + 1.0_dp)**(-1.0_dp/lambda)
      if (.not. lower_tail) p = 1.0_dp - p
   end function paranda_ordaz

   pure elemental real(dp) function daranda_ordaz(eta, lambda, give_log) result(y)
      real(dp), intent(in) :: eta !! Linear predictor for the Aranda-Ordaz density.
      real(dp), intent(in) :: lambda !! Positive Aranda-Ordaz shape parameter.
      logical, intent(in) :: give_log !! Return the log-density when true.
      real(dp) :: ly
      if (lambda <= 0.0_dp) then
         y = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      ly = eta - (1.0_dp + 1.0_dp/lambda)*log(lambda*exp(eta) + 1.0_dp)
      y = merge(ly, exp(ly), give_log)
   end function daranda_ordaz

   pure elemental real(dp) function garanda_ordaz(eta, lambda) result(g)
      real(dp), intent(in) :: eta !! Linear predictor for the derivative of the Aranda-Ordaz density.
      real(dp), intent(in) :: lambda !! Positive Aranda-Ordaz shape parameter.
      real(dp) :: lex, y
      if (lambda <= 0.0_dp) then
         g = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      lex = lambda*exp(eta)
      y = daranda_ordaz(eta, lambda, .false.)
      g = y*(1.0_dp - (1.0_dp + 1.0_dp/lambda)*lex/(1.0_dp + lex))
   end function garanda_ordaz

   pure elemental real(dp) function plgamma_link(eta, lambda, lower_tail) result(p)
      real(dp), intent(in) :: eta !! Linear predictor for the log-gamma link CDF.
      real(dp), intent(in) :: lambda !! Log-gamma shape parameter; zero gives the probit limit.
      logical, intent(in) :: lower_tail !! Return the lower tail when true and upper tail otherwise.
      real(dp) :: q2, v
      if (abs(lambda) <= 1.0e-6_dp) then
         p = normal_cdf(eta)
      else
         q2 = 1.0_dp/(lambda*lambda)
         v = q2*exp(lambda*eta)
         if (lambda > 0.0_dp) then
            p = regularized_gamma_p(q2, v)
         else
            p = 1.0_dp - regularized_gamma_p(q2, v)
         end if
      end if
      if (.not. lower_tail) p = 1.0_dp - p
   end function plgamma_link

   pure elemental real(dp) function dlgamma_link(x, lambda, give_log) result(y)
      real(dp), intent(in) :: x !! Linear predictor for the log-gamma density.
      real(dp), intent(in) :: lambda !! Log-gamma shape parameter; zero gives the normal limit.
      logical, intent(in) :: give_log !! Return the log-density when true.
      real(dp) :: q2, z, ly
      if (abs(lambda) < 1.0e-5_dp) then
         ly = log(normal_pdf(x))
      else
         q2 = 1.0_dp/(lambda*lambda)
         z = lambda*x
         ly = log(abs(lambda)) + q2*log(q2) - log_gamma(q2) + q2*(z - exp(z))
      end if
      y = merge(ly, exp(ly), give_log)
   end function dlgamma_link

   pure elemental real(dp) function glgamma_link(x, lambda) result(g)
      real(dp), intent(in) :: x !! Linear predictor for the derivative of the log-gamma density.
      real(dp), intent(in) :: lambda !! Log-gamma shape parameter; zero gives the normal limit.
      real(dp) :: y
      if (abs(lambda) < 1.0e-5_dp) then
         g = -x*normal_pdf(x)
      else
         y = dlgamma_link(x, lambda, .false.)
         g = y*(1.0_dp - exp(lambda*x))/lambda
      end if
   end function glgamma_link

   pure elemental real(dp) function link_cdf(x, link, lambda, lower_tail) result(p)
      real(dp), intent(in) :: x !! Linear predictor at which the selected link CDF is evaluated.
      integer, intent(in) :: link !! Link identifier from link_logit through link_log_gamma.
      real(dp), intent(in) :: lambda !! Shape parameter used by flexible links and ignored otherwise.
      logical, intent(in) :: lower_tail !! Return the lower tail when true and the upper tail otherwise.
      select case (link)
      case (link_logit)
         if (x >= 0.0_dp) then
            p = 1.0_dp/(1.0_dp + exp(-x))
         else
            p = exp(x)/(1.0_dp + exp(x))
         end if
         if (.not. lower_tail) p = 1.0_dp - p
      case (link_probit)
         p = normal_cdf(x)
         if (.not. lower_tail) p = 1.0_dp - p
      case (link_cloglog)
         p = pgumbel(x, 0.0_dp, 1.0_dp, lower_tail, .false.)
      case (link_loglog)
         p = pgumbel(x, 0.0_dp, 1.0_dp, lower_tail, .true.)
      case (link_cauchit)
         p = 0.5_dp + atan(x)/pi
         if (.not. lower_tail) p = 1.0_dp - p
      case (link_aranda_ordaz)
         p = paranda_ordaz(x, lambda, lower_tail)
      case (link_log_gamma)
         p = plgamma_link(x, lambda, lower_tail)
      case default
         p = ieee_value(0.0_dp, ieee_quiet_nan)
      end select
   end function link_cdf

   pure elemental real(dp) function link_pdf(x, link, lambda) result(y)
      real(dp), intent(in) :: x !! Linear predictor at which the selected link density is evaluated.
      integer, intent(in) :: link !! Link identifier from link_logit through link_log_gamma.
      real(dp), intent(in) :: lambda !! Shape parameter used by flexible links and ignored otherwise.
      real(dp) :: p
      select case (link)
      case (link_logit)
         p = link_cdf(x, link_logit, lambda, .true.)
         y = p*(1.0_dp - p)
      case (link_probit)
         y = normal_pdf(x)
      case (link_cloglog)
         y = dgumbel(x, 0.0_dp, 1.0_dp, .false., .false.)
      case (link_loglog)
         y = dgumbel(x, 0.0_dp, 1.0_dp, .false., .true.)
      case (link_cauchit)
         y = 1.0_dp/(pi*(1.0_dp + x*x))
      case (link_aranda_ordaz)
         y = daranda_ordaz(x, lambda, .false.)
      case (link_log_gamma)
         y = dlgamma_link(x, lambda, .false.)
      case default
         y = ieee_value(0.0_dp, ieee_quiet_nan)
      end select
   end function link_pdf

   pure elemental real(dp) function link_pdf_gradient(x, link, lambda) result(g)
      real(dp), intent(in) :: x !! Linear predictor at which the selected density derivative is evaluated.
      integer, intent(in) :: link !! Link identifier from link_logit through link_log_gamma.
      real(dp), intent(in) :: lambda !! Shape parameter used by flexible links and ignored otherwise.
      select case (link)
      case (link_logit)
         g = glogis(x)
      case (link_probit)
         g = gnorm(x)
      case (link_cloglog)
         g = ggumbel(x, .false.)
      case (link_loglog)
         g = ggumbel(x, .true.)
      case (link_cauchit)
         g = gcauchy(x)
      case (link_aranda_ordaz)
         g = garanda_ordaz(x, lambda)
      case (link_log_gamma)
         g = glgamma_link(x, lambda)
      case default
         g = ieee_value(0.0_dp, ieee_quiet_nan)
      end select
   end function link_pdf_gradient
end module ordinal_links
