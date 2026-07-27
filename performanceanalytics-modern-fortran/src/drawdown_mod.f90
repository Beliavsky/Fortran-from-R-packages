! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module drawdown_mod
  use kinds_mod, only: dp
  use statistics_mod, only: mean_value, sd_value, quantile_type7
  implicit none
  private
  public :: drawdown_series, max_drawdown, drawdown_episode, find_drawdowns
  public :: average_drawdown, average_drawdown_length, average_recovery_length
  public :: drawdown_deviation, pain_index, ulcer_index, conditional_drawdown_at_risk
  public :: drawdown_peak
  type :: drawdown_episode
    integer :: from_index=0
    integer :: trough_index=0
    integer :: to_index=0
    integer :: length=0
    integer :: recovery=0
    real(dp) :: depth=0.0_dp
  end type drawdown_episode
contains
  pure subroutine drawdown_series(r,dd,geometric)
    real(dp),intent(in)::r(:)
    real(dp),intent(out)::dd(:)
    logical,intent(in),optional::geometric
    logical::g
    real(dp)::wealth,peak
    integer::i,n
    n=min(size(r),size(dd)); g=.true.; if(present(geometric)) g=geometric
    wealth=1.0_dp; peak=1.0_dp
    do i=1,n
      if(g) then; wealth=wealth*(1.0_dp+r(i)); else; wealth=wealth+r(i); end if
      peak=max(peak,wealth)
      if(abs(peak)>tiny(1.0_dp)) then; dd(i)=wealth/peak-1.0_dp; else; dd(i)=0.0_dp; end if
    end do
  end subroutine drawdown_series

  pure real(dp) function max_drawdown(r,geometric,invert) result(v)
    real(dp),intent(in)::r(:)
    logical,intent(in),optional::geometric,invert
    logical::inv
    real(dp),allocatable::dd(:)
    allocate(dd(size(r)))
    if(present(geometric)) then; call drawdown_series(r,dd,geometric); else; call drawdown_series(r,dd); end if
    if(size(dd)==0) then; v=0.0_dp; else; v=minval(dd); end if
    inv=.true.; if(present(invert)) inv=invert
    if(inv) v=-v
  end function max_drawdown

  subroutine find_drawdowns(r,episodes,n_episodes,geometric)
    real(dp),intent(in)::r(:)
    type(drawdown_episode),allocatable,intent(out)::episodes(:)
    integer,intent(out)::n_episodes
    logical,intent(in),optional::geometric
    real(dp),allocatable::dd(:)
    type(drawdown_episode),allocatable::tmp(:)
    integer::n,i,start,trough,count_ep
    n=size(r); allocate(dd(n),tmp(max(1,n)))
    if(present(geometric)) then; call drawdown_series(r,dd,geometric); else; call drawdown_series(r,dd); end if
    count_ep=0; i=1
    do while(i<=n)
      if(dd(i)<-epsilon(1.0_dp)) then
        start=max(1,i-1); trough=i
        do while(i<=n)
          if(dd(i)>=-epsilon(1.0_dp)) exit
          if(dd(i)<dd(trough)) trough=i
          i=i+1
        end do
        count_ep=count_ep+1
        tmp(count_ep)%from_index=start
        tmp(count_ep)%trough_index=trough
        if(i<=n) then; tmp(count_ep)%to_index=i; else; tmp(count_ep)%to_index=n; end if
        tmp(count_ep)%length=tmp(count_ep)%to_index-start
        tmp(count_ep)%recovery=tmp(count_ep)%to_index-trough
        tmp(count_ep)%depth=dd(trough)
      else
        i=i+1
      end if
    end do
    n_episodes=count_ep; allocate(episodes(count_ep))
    if(count_ep>0) episodes=tmp(:count_ep)
  end subroutine find_drawdowns

  real(dp) function average_drawdown(r) result(v)
    real(dp),intent(in)::r(:)
    type(drawdown_episode),allocatable::ep(:)
    integer::n,i
    call find_drawdowns(r,ep,n)
    if(n==0) then; v=0.0_dp; else; v=sum([(-ep(i)%depth,i=1,n)])/real(n,dp); end if
  end function average_drawdown

  real(dp) function average_drawdown_length(r) result(v)
    real(dp),intent(in)::r(:)
    type(drawdown_episode),allocatable::ep(:)
    integer::n,i
    call find_drawdowns(r,ep,n)
    if(n==0) then; v=0.0_dp; else; v=real(sum([(ep(i)%length,i=1,n)]),dp)/real(n,dp); end if
  end function average_drawdown_length

  real(dp) function average_recovery_length(r) result(v)
    real(dp),intent(in)::r(:)
    type(drawdown_episode),allocatable::ep(:)
    integer::n,i
    call find_drawdowns(r,ep,n)
    if(n==0) then; v=0.0_dp; else; v=real(sum([(ep(i)%recovery,i=1,n)]),dp)/real(n,dp); end if
  end function average_recovery_length

  real(dp) function drawdown_deviation(r) result(v)
    real(dp),intent(in)::r(:)
    type(drawdown_episode),allocatable::ep(:)
    real(dp),allocatable::x(:)
    integer::n,i
    call find_drawdowns(r,ep,n)
    if(n<=1) then; v=0.0_dp; return; end if
    allocate(x(n)); x=[(-ep(i)%depth,i=1,n)]; v=sd_value(x)
  end function drawdown_deviation

  pure real(dp) function pain_index(r) result(v)
    real(dp),intent(in)::r(:)
    real(dp),allocatable::dd(:)
    allocate(dd(size(r))); call drawdown_series(r,dd)
    if(size(r)==0) then; v=0.0_dp; else; v=-mean_value(dd); end if
  end function pain_index

  pure real(dp) function ulcer_index(r) result(v)
    real(dp),intent(in)::r(:)
    real(dp),allocatable::dd(:)
    allocate(dd(size(r))); call drawdown_series(r,dd)
    if(size(r)==0) then; v=0.0_dp; else; v=sqrt(mean_value(dd*dd)); end if
  end function ulcer_index

  real(dp) function conditional_drawdown_at_risk(r,p) result(v)
    real(dp),intent(in)::r(:),p
    real(dp),allocatable::dd(:),loss(:)
    real(dp)::q
    integer::n
    allocate(dd(size(r))); call drawdown_series(r,dd); loss=-dd
    q=quantile_type7(loss,p); n=count(loss>=q)
    if(n==0) then; v=q; else; v=sum(loss,mask=loss>=q)/real(n,dp); end if
  end function conditional_drawdown_at_risk

  pure integer function drawdown_peak(r) result(idx)
    real(dp),intent(in)::r(:)
    real(dp),allocatable::dd(:)
    integer::loc(1)
    allocate(dd(size(r))); call drawdown_series(r,dd)
    if(size(r)==0) then; idx=0; else; loc=minloc(dd); idx=loc(1); end if
  end function drawdown_peak
end module drawdown_mod
