! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013-2026 the original stochvol and factorstochvol authors
! and the modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2, or (at your option) any later version.
! This file is distributed without any warranty; see LICENSE for details.
module sv_rng
  use sv_kinds, only : dp, pi
  implicit none
  private
  public :: seed_rng, randu, randn, rand_gamma, rand_inv_gamma, rand_student_t, rand_gig_slice
contains
  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i=1,n
      put(i)=modulo(seed + 104729*i + 8191*i*i, huge(1)-1)
      if (put(i)==0) put(i)=i
    end do
    call random_seed(put=put)
  end subroutine seed_rng

  real(dp) function randu() result(x)
    call random_number(x)
    x=max(x,epsilon(1.0_dp))
  end function randu

  real(dp) function randn() result(z)
    real(dp) :: u1,u2
    u1=randu(); u2=randu()
    z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function randn

  recursive real(dp) function rand_gamma(shape, rate) result(x)
    real(dp), intent(in) :: shape, rate
    real(dp) :: d,c,z,u,v
    if (shape<=0.0_dp .or. rate<=0.0_dp) error stop 'rand_gamma: invalid parameters'
    if (shape<1.0_dp) then
      x=rand_gamma(shape+1.0_dp,rate)*randu()**(1.0_dp/shape)
      return
    end if
    d=shape-1.0_dp/3.0_dp
    c=1.0_dp/sqrt(9.0_dp*d)
    do
      z=randn(); v=(1.0_dp+c*z)**3
      if (v<=0.0_dp) cycle
      u=randu()
      if (u<1.0_dp-0.0331_dp*z**4) exit
      if (log(u)<0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
    end do
    x=d*v/rate
  end function rand_gamma

  real(dp) function rand_inv_gamma(shape, scale) result(x)
    real(dp), intent(in) :: shape, scale
    x=1.0_dp/rand_gamma(shape,scale)
  end function rand_inv_gamma

  real(dp) function rand_student_t(nu, standardized) result(x)
    real(dp), intent(in) :: nu
    logical, intent(in), optional :: standardized
    logical :: std
    std=.false.; if (present(standardized)) std=standardized
    x=randn()/sqrt(rand_gamma(0.5_dp*nu,0.5_dp*nu))
    if (std .and. nu>2.0_dp) x=x*sqrt((nu-2.0_dp)/nu)
  end function rand_student_t

  real(dp) function rand_gig_slice(lambda,chi,psi,current) result(x)
    real(dp), intent(in) :: lambda,chi,psi,current
    real(dp) :: z, logy, left, right, proposal
    integer :: k
    if (chi<0.0_dp .or. psi<=0.0_dp) error stop 'rand_gig_slice: invalid parameters'
    z=log(max(current,1.0e-12_dp))
    logy=gig_logz(z,lambda,chi,psi)+log(randu())
    left=z-randu(); right=left+1.0_dp
    do k=1,100
      if (gig_logz(left,lambda,chi,psi)<=logy) exit
      left=left-1.0_dp
    end do
    do k=1,100
      if (gig_logz(right,lambda,chi,psi)<=logy) exit
      right=right+1.0_dp
    end do
    do k=1,10000
      proposal=left+(right-left)*randu()
      if (gig_logz(proposal,lambda,chi,psi)>=logy) then
        x=exp(proposal)
        return
      end if
      if (proposal<z) then
        left=proposal
      else
        right=proposal
      end if
    end do
    x=exp(z)
  end function rand_gig_slice

  pure real(dp) function gig_logz(z,lambda,chi,psi) result(v)
    real(dp), intent(in) :: z,lambda,chi,psi
    v=lambda*z-0.5_dp*(chi*exp(-z)+psi*exp(z))
  end function gig_logz
end module sv_rng
