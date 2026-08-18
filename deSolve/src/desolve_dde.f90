! SPDX-License-Identifier: GPL-2.0-or-later
module desolve_dde
  use desolve_kinds, only : dp
  use desolve_types, only : ode_result
  use desolve_utilities, only : history_buffer
  implicit none
  private
  abstract interface
    subroutine dde_rhs(t,y,history,dydt)
      import dp,history_buffer
      real(dp),intent(in)::t,y(:)
      type(history_buffer),intent(in)::history
      real(dp),intent(out)::dydt(:)
    end subroutine dde_rhs
  end interface
  public :: dde_rhs, dede_rk4
contains
  function dede_rk4(rhs,y0,times,h,initial_history,max_steps) result(sol)
    procedure(dde_rhs)::rhs
    real(dp),intent(in)::y0(:),times(:)
    real(dp),intent(in),optional::h
    type(history_buffer),intent(in),optional::initial_history
    integer,intent(in),optional::max_steps
    type(ode_result)::sol
    type(history_buffer)::hist
    real(dp),allocatable::y(:),yn(:),k1(:),k2(:),k3(:),k4(:),yt(:)
    real(dp)::t,target,hs,hbase,dir
    integer::i,n,nst,maxst
    n=size(y0);allocate(sol%t(size(times)),sol%y(n,size(times)));sol%t=times;sol%y=0.0_dp;sol%y(:,1)=y0
    if(.not.valid_times(times))then;sol%status=-100;sol%message='times must be strictly monotone';return;end if
    if(present(initial_history))then;hist=initial_history;else;call hist%init(n,max(1000,100*size(times)));end if
    allocate(y(n),yn(n),k1(n),k2(n),k3(n),k4(n),yt(n));y=y0;t=times(1)
    if(hist%count==0)then
      k1=0.0_dp;call hist%append(t,y,k1);call rhs(t,y,hist,k1);hist%dy(:,hist%count)=k1
    else if(abs(hist%time(hist%count)-t)>100.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(t)))then
      call rhs(t,y,hist,k1);call hist%append(t,y,k1)
    end if
    dir=merge(1.0_dp,-1.0_dp,times(size(times))>times(1));hbase=dir*abs(times(size(times))-times(1))/1000.0_dp
    if(present(h))hbase=dir*abs(h);if(abs(hbase)<=tiny(1.0_dp))hbase=dir*1e-4_dp
    maxst=100000;if(present(max_steps))maxst=max_steps;nst=0
    do i=2,size(times)
      target=times(i)
      do while(dir*(target-t)>100.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(target)))
        if(nst>=maxst)then;sol%status=-1;sol%message='maximum DDE steps exceeded';return;end if
        hs=hbase;if(dir*(t+hs-target)>0.0_dp)hs=target-t
        call rhs(t,y,hist,k1);yt=y+0.5_dp*hs*k1;call rhs(t+0.5_dp*hs,yt,hist,k2)
        yt=y+0.5_dp*hs*k2;call rhs(t+0.5_dp*hs,yt,hist,k3);yt=y+hs*k3;call rhs(t+hs,yt,hist,k4)
        yn=y+hs*(k1+2.0_dp*k2+2.0_dp*k3+k4)/6.0_dp;t=t+hs;y=yn;call rhs(t,y,hist,k1);call hist%append(t,y,k1)
        nst=nst+1
      end do
      sol%y(:,i)=y
    end do
    sol%status=0;sol%message='success';sol%stats%n_steps=nst;sol%stats%n_rhs=5*nst+1;sol%stats%step_last=hbase
  end function dede_rk4
  pure logical function valid_times(t)result(ok)
    real(dp),intent(in)::t(:);integer::i;ok=size(t)>=1;if(size(t)<2)return
    if(t(2)>t(1))then;do i=2,size(t);if(t(i)<=t(i-1))then;ok=.false.;return;end if;end do
    else if(t(2)<t(1))then;do i=2,size(t);if(t(i)>=t(i-1))then;ok=.false.;return;end if;end do
    else;ok=.false.;end if
  end function valid_times
end module desolve_dde
