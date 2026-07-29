! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2018 Marius Hofert, Erik Hintz and Christiane Lemieux
module nvmix_gamma_mix
  use nvmix_kinds, only : dp, i8
  use nvmix_types
  use nvmix_special, only : chi_square_pdf,chi_square_cdf,chi_square_quantile,f_pdf,f_cdf
  use nvmix_random, only : halton,seed_random,random_gamma
  use nvmix_mixing, only : mixing_quantile,mixing_random
  implicit none
  private
  public :: dgammamix,pgammamix,qgammamix,rgammamix
contains
  real(dp) function dgammamix(x,d,family,parameter,control,log_density) result(v)
    real(dp), intent(in) :: x,parameter
    integer, intent(in) :: d,family
    type(integration_control), intent(in), optional :: control
    logical, intent(in), optional :: log_density
    type(integration_control) :: ctrl
    logical :: lg
    real(dp) :: w,s
    integer :: i,n
    ctrl=integration_control(); if(present(control))ctrl=control
    lg=.false.; if(present(log_density))lg=log_density
    if(x<=0.0_dp .or. d<1)then; v=0.0_dp; if(lg)v=-huge(1.0_dp); return; end if
    select case(family)
    case(mix_constant); v=chi_square_pdf(x,real(d,dp))
    case(mix_inverse_gamma); v=f_pdf(x/real(d,dp),real(d,dp),parameter)/real(d,dp)
    case default
      n=max(256,ctrl%samples); s=0.0_dp
      do i=1,n
        w=mixing_quantile(halton(i,2),family,parameter)
        s=s+chi_square_pdf(x/w,real(d,dp))/w
      end do
      v=s/real(n,dp)
    end select
    if(lg)v=log(max(v,tiny(1.0_dp)))
  end function
  real(dp) function pgammamix(x,d,family,parameter,control) result(v)
    real(dp), intent(in) :: x,parameter
    integer, intent(in) :: d,family
    type(integration_control), intent(in), optional :: control
    type(integration_control) :: ctrl
    real(dp) :: w,s
    integer :: i,n
    ctrl=integration_control(); if(present(control))ctrl=control
    if(x<=0.0_dp)then; v=0.0_dp; return; end if
    select case(family)
    case(mix_constant); v=chi_square_cdf(x,real(d,dp))
    case(mix_inverse_gamma); v=f_cdf(x/real(d,dp),real(d,dp),parameter)
    case default
      n=max(256,ctrl%samples); s=0.0_dp
      do i=1,n
        w=mixing_quantile(halton(i,2),family,parameter)
        s=s+chi_square_cdf(x/w,real(d,dp))
      end do
      v=s/real(n,dp)
    end select
    v=min(1.0_dp,max(0.0_dp,v))
  end function
  real(dp) function qgammamix(p,d,family,parameter,control) result(v)
    real(dp), intent(in) :: p,parameter
    integer, intent(in) :: d,family
    type(integration_control), intent(in), optional :: control
    type(integration_control) :: ctrl
    real(dp) :: lo,hi,mid
    integer :: i
    ctrl=integration_control(); if(present(control))ctrl=control
    if(p<=0.0_dp)then; v=0.0_dp; return; end if
    if(p>=1.0_dp)then; v=huge(1.0_dp); return; end if
    lo=0.0_dp; hi=max(1.0_dp,real(d,dp))
    do while(pgammamix(hi,d,family,parameter,ctrl)<p); hi=2.0_dp*hi; end do
    do i=1,100
      mid=0.5_dp*(lo+hi)
      if(pgammamix(mid,d,family,parameter,ctrl)<p)then; lo=mid; else; hi=mid; end if
    end do
    v=0.5_dp*(lo+hi)
  end function
  function rgammamix(n,d,family,parameter,seed) result(values)
    integer, intent(in) :: n,d,family
    real(dp), intent(in) :: parameter
    integer(i8), intent(in), optional :: seed
    real(dp), allocatable :: values(:)
    integer :: i
    if(present(seed))call seed_random(seed)
    allocate(values(max(0,n)))
    do i=1,n; values(i)=mixing_random(family,parameter)*random_gamma(0.5_dp*real(d,dp),2.0_dp); end do
  end function
end module nvmix_gamma_mix
