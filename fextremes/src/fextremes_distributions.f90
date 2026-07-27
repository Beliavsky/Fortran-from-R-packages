! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software: you can redistribute it
! and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 2 of the License,
! or (at your option) any later version.
module fextremes_distributions
  use fextremes_kinds, only: dp, pi, euler_gamma
  use fextremes_rng, only: rng_state, uniform_rng, exponential_rng
  use fextremes_stats, only: nan_value
  implicit none
  private
  public :: gev_pdf, gev_logpdf, gev_cdf, gev_quantile, gev_random, gev_sample, gev_moments
  public :: gpd_pdf, gpd_logpdf, gpd_cdf, gpd_quantile, gpd_random, gpd_sample, gpd_moments
contains
  pure real(dp) function gev_logpdf(x, xi, mu, beta) result(v)
    real(dp), intent(in) :: x, xi, mu, beta
    real(dp) :: z, t
    if (beta <= 0.0_dp) then; v=-huge(1.0_dp); return; end if
    z=(x-mu)/beta
    if(abs(xi)<1.0e-10_dp) then
      v=-log(beta)-z-exp(-z)
    else
      t=1.0_dp+xi*z
      if(t<=0.0_dp) then; v=-huge(1.0_dp)
      else; v=-log(beta)-t**(-1.0_dp/xi)-(1.0_dp/xi+1.0_dp)*log(t); end if
    end if
  end function gev_logpdf

  pure real(dp) function gev_pdf(x, xi, mu, beta) result(v)
    real(dp), intent(in) :: x, xi, mu, beta
    real(dp) :: lv
    lv=gev_logpdf(x,xi,mu,beta)
    if(lv < log(tiny(1.0_dp))) then; v=0.0_dp; else; v=exp(lv); end if
  end function gev_pdf

  pure real(dp) function gev_cdf(x, xi, mu, beta) result(p)
    real(dp), intent(in) :: x, xi, mu, beta
    real(dp) :: z,t
    if(beta<=0.0_dp) then; p=nan_value(); return; end if
    z=(x-mu)/beta
    if(abs(xi)<1.0e-10_dp) then
      p=exp(-exp(-z))
    else
      t=1.0_dp+xi*z
      if(t<=0.0_dp) then
        if(xi>0.0_dp) then; p=0.0_dp; else; p=1.0_dp; end if
      else
        p=exp(-t**(-1.0_dp/xi))
      end if
    end if
  end function gev_cdf

  pure real(dp) function gev_quantile(p, xi, mu, beta) result(q)
    real(dp), intent(in) :: p, xi, mu, beta
    if(beta<=0.0_dp .or. p<0.0_dp .or. p>1.0_dp) then; q=nan_value(); return; end if
    if(p<=0.0_dp) then
      if(xi>0.0_dp) then; q=mu-beta/xi; else; q=-huge(1.0_dp); end if
    else if(p>=1.0_dp) then
      if(xi<0.0_dp) then; q=mu-beta/xi; else; q=huge(1.0_dp); end if
    else if(abs(xi)<1.0e-10_dp) then
      q=mu-beta*log(-log(p))
    else
      q=mu+beta*((-log(p))**(-xi)-1.0_dp)/xi
    end if
  end function gev_quantile

  real(dp) function gev_random(rng, xi, mu, beta) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: xi,mu,beta
    x=gev_quantile(uniform_rng(rng),xi,mu,beta)
  end function gev_random

  subroutine gev_sample(rng,xi,mu,beta,x)
    type(rng_state),intent(inout)::rng
    real(dp),intent(in)::xi,mu,beta
    real(dp),intent(out)::x(:)
    integer::i
    do i=1,size(x); x(i)=gev_random(rng,xi,mu,beta); end do
  end subroutine gev_sample

  subroutine gev_moments(xi,mu,beta,mean,var)
    real(dp), intent(in) :: xi,mu,beta
    real(dp), intent(out) :: mean,var
    real(dp) :: g1,g2
    mean=nan_value(); var=nan_value()
    if(beta<=0.0_dp) return
    if(abs(xi)<1.0e-10_dp) then
      mean=mu+beta*euler_gamma; var=(beta*pi)**2/6.0_dp; return
    end if
    if(xi<1.0_dp) then
      g1=gamma(1.0_dp-xi); mean=mu+beta*(g1-1.0_dp)/xi
    end if
    if(xi<0.5_dp) then
      g1=gamma(1.0_dp-xi); g2=gamma(1.0_dp-2.0_dp*xi)
      var=(beta/xi)**2*(g2-g1*g1)
    end if
  end subroutine gev_moments

  pure real(dp) function gpd_logpdf(x, xi, mu, beta) result(v)
    real(dp), intent(in) :: x,xi,mu,beta
    real(dp) :: z,t
    if(beta<=0.0_dp .or. x<mu) then; v=-huge(1.0_dp); return; end if
    z=(x-mu)/beta
    if(abs(xi)<1.0e-10_dp) then
      v=-log(beta)-z
    else
      t=1.0_dp+xi*z
      if(t<=0.0_dp) then; v=-huge(1.0_dp); else; v=-log(beta)-(1.0_dp/xi+1.0_dp)*log(t); end if
    end if
  end function gpd_logpdf

  pure real(dp) function gpd_pdf(x, xi, mu, beta) result(v)
    real(dp), intent(in) :: x,xi,mu,beta
    real(dp) :: lv
    lv=gpd_logpdf(x,xi,mu,beta)
    if(lv<log(tiny(1.0_dp))) then; v=0.0_dp; else; v=exp(lv); end if
  end function gpd_pdf

  pure real(dp) function gpd_cdf(x, xi, mu, beta) result(p)
    real(dp), intent(in) :: x,xi,mu,beta
    real(dp) :: z,t
    if(beta<=0.0_dp) then; p=nan_value(); return; end if
    if(x<=mu) then; p=0.0_dp; return; end if
    z=(x-mu)/beta
    if(abs(xi)<1.0e-10_dp) then
      p=1.0_dp-exp(-z)
    else
      t=1.0_dp+xi*z
      if(t<=0.0_dp) then; p=1.0_dp; else; p=1.0_dp-t**(-1.0_dp/xi); end if
    end if
  end function gpd_cdf

  pure real(dp) function gpd_quantile(p, xi, mu, beta) result(q)
    real(dp), intent(in) :: p,xi,mu,beta
    if(beta<=0.0_dp .or. p<0.0_dp .or. p>1.0_dp) then; q=nan_value(); return; end if
    if(p<=0.0_dp) then; q=mu; return; end if
    if(p>=1.0_dp) then
      if(xi<0.0_dp) then; q=mu-beta/xi; else; q=huge(1.0_dp); end if
    else if(abs(xi)<1.0e-10_dp) then
      q=mu-beta*log(1.0_dp-p)
    else
      q=mu+beta*((1.0_dp-p)**(-xi)-1.0_dp)/xi
    end if
  end function gpd_quantile

  real(dp) function gpd_random(rng, xi, mu, beta) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: xi,mu,beta
    x=gpd_quantile(uniform_rng(rng),xi,mu,beta)
  end function gpd_random

  subroutine gpd_sample(rng,xi,mu,beta,x)
    type(rng_state),intent(inout)::rng
    real(dp),intent(in)::xi,mu,beta
    real(dp),intent(out)::x(:)
    integer::i
    do i=1,size(x); x(i)=gpd_random(rng,xi,mu,beta); end do
  end subroutine gpd_sample

  subroutine gpd_moments(xi,mu,beta,mean,var)
    real(dp), intent(in) :: xi,mu,beta
    real(dp), intent(out) :: mean,var
    mean=nan_value(); var=nan_value()
    if(beta<=0.0_dp) return
    if(xi<1.0_dp) mean=mu+beta/(1.0_dp-xi)
    if(xi<0.5_dp) var=beta*beta/((1.0_dp-xi)**2*(1.0_dp-2.0_dp*xi))
  end subroutine gpd_moments
end module fextremes_distributions
