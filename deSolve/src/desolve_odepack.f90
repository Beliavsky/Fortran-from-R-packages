! SPDX-License-Identifier: GPL-2.0-or-later
module desolve_odepack
  use desolve_kinds, only : dp
  use desolve_types, only : ode_rhs, ode_result
  use desolve_callback_bridge, only : current_rhs, odepack_rhs_bridge, &
      odepack_jac_dummy, vode_rhs_bridge, vode_jac_dummy
  implicit none
  private
  public :: lsoda, lsode, vode

  interface
    subroutine dlsoda(f, neq, y, t, tout, itol, rtol, atol, itask, istate, &
        iopt, rwork, lrw, iwork, liw, jac, jt, rpar, ipar)
      import dp
      external :: f, jac
      integer :: neq(*), itol, itask, istate, iopt, lrw, iwork(*), liw, jt, ipar(*)
      real(dp) :: y(*), t, tout, rtol(*), atol(*), rwork(*), rpar(*)
    end subroutine dlsoda

    subroutine dlsode(f, neq, y, t, tout, itol, rtol, atol, itask, istate, &
        iopt, rwork, lrw, iwork, liw, jac, mf, rpar, ipar)
      import dp
      external :: f, jac
      integer :: neq(*), itol, itask, istate, iopt, lrw, iwork(*), liw, mf, ipar(*)
      real(dp) :: y(*), t, tout, rtol(*), atol(*), rwork(*), rpar(*)
    end subroutine dlsode

    subroutine dvode(f, neq, y, t, tout, itol, rtol, atol, itask, istate, &
        iopt, rwork, lrw, iwork, liw, jac, mf, rpar, ipar)
      import dp
      external :: f, jac
      integer :: neq, itol, itask, istate, iopt, lrw, iwork(*), liw, mf, ipar(*)
      real(dp) :: y(*), t, tout, rtol(*), atol(*), rwork(*), rpar(*)
    end subroutine dvode
  end interface

contains

  function lsoda(rhs, y0, times, rtol, atol, max_steps) result(sol)
    procedure(ode_rhs) :: rhs
    real(dp), intent(in) :: y0(:), times(:)
    real(dp), intent(in), optional :: rtol, atol
    integer, intent(in), optional :: max_steps
    type(ode_result) :: sol
    integer :: n, nt, neq(1), itol, itask, istate, iopt, lrw, liw, jt, i
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: rwork(:), y(:)
    real(dp) :: rt(1), at(1), t, rpar(1)
    integer :: ipar(1)

    n=size(y0); nt=size(times)
    call init_result(sol,n,nt,times,y0)
    if (.not. valid_times(times)) then
      sol%status=-100; sol%message='times must be strictly monotone'; return
    end if
    if (n < 1 .or. nt < 1) then
      sol%status=-101; sol%message='empty state or time vector'; return
    end if

    neq(1)=n; itol=1; itask=1; istate=1; iopt=0; jt=2
    lrw=max(100, 22+n*max(16,n+9)+20*n)
    liw=max(50,20+n)
    allocate(rwork(lrw),iwork(liw),y(n)); rwork=0.0_dp; iwork=0
    y=y0; t=times(1); rt(1)=1.0e-8_dp; at(1)=1.0e-10_dp
    if (present(rtol)) rt(1)=rtol
    if (present(atol)) at(1)=atol
    if (present(max_steps)) then
      iopt=1; iwork(6)=max_steps
    end if
    rpar=0.0_dp; ipar=0
    current_rhs => rhs
    do i=2,nt
      call dlsoda(odepack_rhs_bridge,neq,y,t,times(i),itol,rt,at,itask,istate,iopt, &
          rwork,lrw,iwork,liw,odepack_jac_dummy,jt,rpar,ipar)
      sol%y(:,i)=y
      if (istate < 0) exit
    end do
    nullify(current_rhs)
    if (istate >= 0 .and. i <= nt) then
      if (i < nt) sol%y(:,i+1:nt)=spread(y,2,nt-i)
    end if
    sol%status=istate
    sol%message=odepack_message(istate)
    call fill_stats(sol,iwork,rwork)
  end function lsoda

  function lsode(rhs, y0, times, rtol, atol, mf, max_steps) result(sol)
    procedure(ode_rhs) :: rhs
    real(dp), intent(in) :: y0(:), times(:)
    real(dp), intent(in), optional :: rtol, atol
    integer, intent(in), optional :: mf, max_steps
    type(ode_result) :: sol
    integer :: n,nt,neq(1),itol,itask,istate,iopt,lrw,liw,mf0,i
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: rwork(:),y(:)
    real(dp) :: rt(1),at(1),t,rpar(1)
    integer :: ipar(1)
    n=size(y0); nt=size(times)
    call init_result(sol,n,nt,times,y0)
    if (.not.valid_times(times)) then
      sol%status=-100; sol%message='times must be strictly monotone'; return
    end if
    neq(1)=n; itol=1; itask=1; istate=1; iopt=0; mf0=22
    if (present(mf)) mf0=mf
    lrw=max(100, 100+30*n+n*n); liw=max(50,30+n)
    allocate(rwork(lrw),iwork(liw),y(n)); rwork=0.0_dp; iwork=0; y=y0
    rt(1)=1e-8_dp; at(1)=1e-10_dp
    if(present(rtol)) rt(1)=rtol
    if(present(atol)) at(1)=atol
    if(present(max_steps)) then; iopt=1; iwork(6)=max_steps; end if
    t=times(1); rpar=0.0_dp; ipar=0; current_rhs=>rhs
    do i=2,nt
      call dlsode(odepack_rhs_bridge,neq,y,t,times(i),itol,rt,at,itask,istate,iopt, &
          rwork,lrw,iwork,liw,odepack_jac_dummy,mf0,rpar,ipar)
      sol%y(:,i)=y
      if(istate<0) exit
    end do
    nullify(current_rhs)
    sol%status=istate; sol%message=odepack_message(istate)
    call fill_stats(sol,iwork,rwork)
  end function lsode

  function vode(rhs, y0, times, rtol, atol, mf, max_steps) result(sol)
    procedure(ode_rhs) :: rhs
    real(dp), intent(in) :: y0(:),times(:)
    real(dp), intent(in), optional :: rtol,atol
    integer, intent(in), optional :: mf,max_steps
    type(ode_result) :: sol
    integer :: n,nt,itol,itask,istate,iopt,lrw,liw,mf0,i
    integer, allocatable :: iwork(:)
    real(dp), allocatable :: rwork(:),y(:)
    real(dp)::rt(1),at(1),t,rpar(1)
    integer::ipar(1)
    n=size(y0);nt=size(times);call init_result(sol,n,nt,times,y0)
    if(.not.valid_times(times))then;sol%status=-100;sol%message='times must be strictly monotone';return;end if
    itol=1;itask=1;istate=1;iopt=0;mf0=22;if(present(mf))mf0=mf
    lrw=max(100,100+30*n+n*n);liw=max(50,30+n)
    allocate(rwork(lrw),iwork(liw),y(n));rwork=0.0_dp;iwork=0;y=y0
    rt(1)=1e-8_dp;at(1)=1e-10_dp;if(present(rtol))rt(1)=rtol;if(present(atol))at(1)=atol
    if(present(max_steps))then;iopt=1;iwork(6)=max_steps;end if
    t=times(1);rpar=0.0_dp;ipar=0;current_rhs=>rhs
    do i=2,nt
      call dvode(vode_rhs_bridge,n,y,t,times(i),itol,rt,at,itask,istate,iopt,rwork,lrw, &
          iwork,liw,vode_jac_dummy,mf0,rpar,ipar)
      sol%y(:,i)=y;if(istate<0)exit
    end do
    nullify(current_rhs);sol%status=istate;sol%message=odepack_message(istate)
    call fill_stats(sol,iwork,rwork)
  end function vode

  subroutine init_result(sol,n,nt,times,y0)
    type(ode_result),intent(out)::sol
    integer,intent(in)::n,nt
    real(dp),intent(in)::times(:),y0(:)
    allocate(sol%t(nt),sol%y(n,nt));sol%t=times;sol%y=0.0_dp
    if(nt>0.and.n>0)sol%y(:,1)=y0
    sol%status=0;sol%message='not run'
  end subroutine init_result

  pure logical function valid_times(t) result(ok)
    real(dp),intent(in)::t(:)
    integer::i
    ok=size(t)>=1;if(size(t)<2)return
    if(t(2)>t(1))then
      do i=2,size(t);if(t(i)<=t(i-1))then;ok=.false.;return;end if;end do
    else if(t(2)<t(1))then
      do i=2,size(t);if(t(i)>=t(i-1))then;ok=.false.;return;end if;end do
    else
      ok=.false.
    end if
  end function valid_times

  subroutine fill_stats(sol,iw,rw)
    type(ode_result),intent(inout)::sol
    integer,intent(in)::iw(:);real(dp),intent(in)::rw(:)
    if(size(iw)>=20)then
      sol%stats%n_steps=iw(11);sol%stats%n_rhs=iw(12);sol%stats%n_jac=iw(13)
      sol%stats%order_last=iw(14);sol%stats%method_last=iw(19)
    end if
    if(size(rw)>=15)then
      sol%stats%step_last=rw(11);sol%stats%step_next=rw(12);sol%stats%time_last_switch=rw(15)
    end if
  end subroutine fill_stats

  pure function odepack_message(istate) result(msg)
    integer,intent(in)::istate;character(len=:),allocatable::msg
    select case(istate)
    case(2);msg='success'
    case(-1);msg='excess work on this call'
    case(-2);msg='too much accuracy requested'
    case(-3);msg='illegal input detected'
    case(-4);msg='repeated error-test failures'
    case(-5);msg='repeated convergence failures'
    case(-6);msg='error weight became zero'
    case(-7);msg='workspace insufficient'
    case default;msg='solver status '//itoa(istate)
    end select
  end function odepack_message

  pure function itoa(i) result(s)
    integer,intent(in)::i;character(len=:),allocatable::s;character(len=32)::buf
    write(buf,'(i0)')i;s=trim(buf)
  end function itoa

end module desolve_odepack
