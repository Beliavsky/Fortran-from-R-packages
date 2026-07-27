! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software: you can redistribute it
! and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 2 of the License,
! or (at your option) any later version.
module fextremes_extremal_index
  use fextremes_kinds, only: dp
  use fextremes_rng, only: rng_state, uniform_rng, exponential_rng
  use fextremes_stats, only: quantile_type1, nan_value, order_real
  implicit none
  private
  public :: theta_result, theta_simulate, block_theta, cluster_theta, run_theta, ferro_seg_theta

  type :: theta_result
    real(dp), allocatable :: probabilities(:), thresholds(:), theta(:)
    integer, allocatable :: exceedances(:), clusters(:), run_lengths(:), messages(:)
  end type theta_result
contains
  subroutine theta_simulate(rng, model, n, theta, x)
    type(rng_state), intent(inout) :: rng
    character(len=*), intent(in) :: model
    integer, intent(in) :: n
    real(dp), intent(in) :: theta
    real(dp), intent(out) :: x(n)
    real(dp), allocatable :: eps(:)
    real(dp) :: th
    integer :: i
    if (n <= 0) return
    if (trim(model) == 'pair') then
      th = 0.5_dp
      allocate(eps(n+1)); eps=0.0_dp
      do i=1,n+1; eps(i)=exponential_rng(rng); end do
      do i=1,n; x(i)=max(eps(i),eps(i+1)); end do
    else
      th=max(min(theta,1.0_dp),epsilon(1.0_dp)); allocate(eps(n)); eps=0.0_dp
      do i=1,n; eps(i)=1.0_dp/(-log(uniform_rng(rng))); end do
      x(1)=th*eps(1)
      do i=2,n; x(i)=max((1.0_dp-th)*x(i-1),th*eps(i)); end do
    end if
  end subroutine theta_simulate

  subroutine block_theta(x, block, probabilities, result)
    real(dp),intent(in)::x(:),probabilities(:)
    integer,intent(in)::block
    type(theta_result),intent(out)::result
    integer::k,n,j,b,lo,hi,nexc,ncl
    real(dp)::u
    k=size(x)/block; n=k*block
    call allocate_result(result,size(probabilities),with_clusters=.true.)
    result%probabilities=probabilities
    do j=1,size(probabilities)
      u=quantile_type1(x,probabilities(j)); nexc=count(x(1:n)>u); ncl=0
      do b=1,k; lo=(b-1)*block+1; hi=b*block; if(maxval(x(lo:hi))>u) ncl=ncl+1; end do
      result%thresholds(j)=u; result%exceedances(j)=nexc; result%clusters(j)=ncl
      if(nexc>0 .and. ncl<k) then
        result%theta(j)=(real(k,dp)/real(n,dp))*log(1.0_dp-real(ncl,dp)/real(k,dp))/ &
          log(1.0_dp-real(nexc,dp)/real(n,dp))
      else; result%theta(j)=nan_value(); end if
    end do
  end subroutine block_theta

  subroutine cluster_theta(x, block, probabilities, result)
    real(dp),intent(in)::x(:),probabilities(:)
    integer,intent(in)::block
    type(theta_result),intent(out)::result
    integer::k,n,j,b,lo,hi,nexc,ncl
    real(dp)::u
    k=size(x)/block; n=k*block
    call allocate_result(result,size(probabilities),with_clusters=.true.); result%probabilities=probabilities
    do j=1,size(probabilities)
      u=quantile_type1(x,probabilities(j)); nexc=count(x(1:n)>u); ncl=0
      do b=1,k; lo=(b-1)*block+1; hi=b*block; if(maxval(x(lo:hi))>u) ncl=ncl+1; end do
      result%thresholds(j)=u; result%exceedances(j)=nexc; result%clusters(j)=ncl
      if(nexc>0) then; result%theta(j)=real(ncl,dp)/real(nexc,dp); else; result%theta(j)=nan_value(); end if
    end do
  end subroutine cluster_theta

  subroutine run_theta(x, run_length, probabilities, result)
    real(dp),intent(in)::x(:),probabilities(:)
    integer,intent(in)::run_length
    type(theta_result),intent(out)::result
    integer::j,i,nexc,nsep,last
    real(dp)::u
    call allocate_result(result,size(probabilities)); result%probabilities=probabilities
    do j=1,size(probabilities)
      u=quantile_type1(x,probabilities(j)); nexc=count(x>u); nsep=0; last=0
      do i=1,size(x)
        if(x(i)>u) then
          if(last>0 .and. i-last>run_length) nsep=nsep+1
          last=i
        end if
      end do
      result%thresholds(j)=u; result%exceedances(j)=nexc
      if(nexc>0) then; result%theta(j)=real(nsep,dp)/real(nexc,dp); else; result%theta(j)=nan_value(); end if
    end do
  end subroutine run_theta

  subroutine ferro_seg_theta(x, probabilities, result)
    real(dp),intent(in)::x(:),probabilities(:)
    type(theta_result),intent(out)::result
    integer,allocatable::pos(:),tt(:),ord(:)
    integer::j,i,nexc,k,tval,msg
    real(dp)::u,theta,s1,s2
    call allocate_result(result,size(probabilities),with_clusters=.true.,with_runs=.true.,with_messages=.true.)
    result%probabilities=probabilities
    do j=1,size(probabilities)
      theta=0.0_dp
      u=quantile_type1(x,probabilities(j)); nexc=count(x>u); result%thresholds(j)=u; result%exceedances(j)=nexc
      if(nexc<2) then; result%theta(j)=nan_value(); cycle; end if
      allocate(pos(nexc),tt(nexc-1),ord(nexc-1)); k=0
      do i=1,size(x); if(x(i)>u) then; k=k+1; pos(k)=i; end if; end do
      tt=pos(2:)-pos(:nexc-1); msg=0
      if(any(tt>2)) then
        s1=sum(real(tt-1,dp)); s2=sum(real((tt-1)*(tt-2),dp))
        if(s2>0.0_dp) theta=2.0_dp*s1*s1/(real(nexc-1,dp)*s2)
        msg=100
      else
        s1=sum(real(tt,dp)); s2=sum(real(tt*tt,dp))
        if(s2>0.0_dp) theta=2.0_dp*s1*s1/(real(nexc-1,dp)*s2)
        msg=1
      end if
      if(theta>1.0_dp) then; theta=1.0_dp; msg=msg*10; end if
      theta=max(0.0_dp,theta); k=max(1,ceiling(theta*real(nexc,dp)))
      call order_real(real(tt,dp),ord,ascending=.false.); tval=tt(ord(min(k,size(tt))))
      result%theta(j)=theta; result%clusters(j)=k; result%run_lengths(j)=tval; result%messages(j)=msg
      deallocate(pos,tt,ord)
    end do
  end subroutine ferro_seg_theta

  subroutine allocate_result(r,n,with_clusters,with_runs,with_messages)
    type(theta_result),intent(out)::r
    integer,intent(in)::n
    logical,intent(in),optional::with_clusters,with_runs,with_messages
    logical::a,b,c
    a=.false.;b=.false.;c=.false.
    if(present(with_clusters)) a=with_clusters
    if(present(with_runs)) b=with_runs
    if(present(with_messages)) c=with_messages
    allocate(r%probabilities(n),r%thresholds(n),r%theta(n),r%exceedances(n))
    if(a) allocate(r%clusters(n)); if(b) allocate(r%run_lengths(n)); if(c) allocate(r%messages(n))
  end subroutine allocate_result
end module fextremes_extremal_index
