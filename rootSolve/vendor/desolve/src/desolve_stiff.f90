! SPDX-License-Identifier: GPL-2.0-or-later
module desolve_stiff
  use desolve_kinds, only : dp
  use desolve_types, only : ode_rhs, complex_ode_rhs, dae_residual, ode_result, complex_ode_result
  use desolve_callback_bridge, only : current_rhs, current_complex_rhs, current_residual, current_n, &
      radau_rhs_bridge, radau_jac_dummy, radau_mass_identity, radau_solout_dummy, &
      daspk_res_bridge, daspk_jac_dummy, daspk_psol_dummy, zvode_rhs_bridge, zvode_jac_dummy
  implicit none
  private
  public :: radau, daspk, zvode

  interface
    subroutine radau5(n,fcn,x,y,xend,h,rtol,atol,itol,jac,ijac,mljac,mujac, &
        mas,imas,mlmas,mumas,solout,iout,work,lwork,iwork,liwork,rpar,ipar,idid)
      import dp
      external :: fcn,jac,mas,solout
      integer :: n,itol,ijac,mljac,mujac,imas,mlmas,mumas,iout,lwork,iwork(*),liwork,ipar(*),idid
      real(dp) :: x,y(*),xend,h,rtol(*),atol(*),work(*),rpar(*)
    end subroutine radau5

    subroutine ddaspk(res,neq,t,y,yprime,tout,info,rtol,atol,idid,rwork,lrw, &
        iwork,liw,rpar,ipar,jac,psol)
      import dp
      external :: res,jac,psol
      integer :: neq,info(*),idid,lrw,iwork(*),liw,ipar(*)
      real(dp) :: t,y(*),yprime(*),tout,rtol(*),atol(*),rwork(*),rpar(*)
    end subroutine ddaspk

  end interface

contains

  function radau(rhs,y0,times,rtol,atol,hinit) result(sol)
    procedure(ode_rhs) :: rhs
    real(dp),intent(in)::y0(:),times(:)
    real(dp),intent(in),optional::rtol,atol,hinit
    type(ode_result)::sol
    integer::n,nt,i,itol,ijac,mljac,mujac,imas,mlmas,mumas,iout,lwork,liwork,idid
    integer,allocatable::iwork(:)
    real(dp),allocatable::work(:),y(:)
    real(dp)::x,h,rt(1),at(1),rpar(1)
    integer::ipar(1)
    n=size(y0);nt=size(times);call init_real(sol,n,nt,times,y0)
    if(.not.valid_times(times))then;sol%status=-100;sol%message='times must be strictly monotone';return;end if
    itol=0;ijac=0;mljac=n;mujac=0;imas=0;mlmas=n;mumas=0;iout=0
    lwork=max(100,4*n*n+20*n+100);liwork=max(50,3*n+50)
    allocate(work(lwork),iwork(liwork),y(n));work=0.0_dp;iwork=0;y=y0
    rt=1e-8_dp;at=1e-10_dp;if(present(rtol))rt=rtol;if(present(atol))at=atol
    idid=1;x=times(1);h=1e-4_dp*sign(1.0_dp,times(nt)-times(1));if(present(hinit))h=sign(abs(hinit),h)
    rpar=0.0_dp;ipar=0;current_rhs=>rhs
    do i=2,nt
      call radau5(n,radau_rhs_bridge,x,y,times(i),h,rt,at,itol,radau_jac_dummy,ijac,mljac,mujac, &
          radau_mass_identity,imas,mlmas,mumas,radau_solout_dummy,iout,work,lwork,iwork,liwork, &
          rpar,ipar,idid)
      sol%y(:,i)=y
      if(idid<0)exit
    end do
    nullify(current_rhs);sol%status=idid
    if(idid>0)then;sol%message='success';else;sol%message='RADAU5 failed with IDID='//itoa(idid);end if
    if(liwork>=14)then;sol%stats%n_rhs=iwork(14);sol%stats%n_jac=iwork(15);sol%stats%n_steps=iwork(16);end if
    sol%stats%step_next=h
  end function radau

  function daspk(residual,y0,yprime0,times,rtol,atol) result(sol)
    procedure(dae_residual)::residual
    real(dp),intent(in)::y0(:),yprime0(:),times(:)
    real(dp),intent(in),optional::rtol,atol
    type(ode_result)::sol
    integer::n,nt,i,idid,lrw,liw
    integer::info(25),ipar(1)
    integer,allocatable::iwork(:)
    real(dp),allocatable::rwork(:),y(:),yp(:)
    real(dp)::t,rt(1),at(1),rpar(1)
    n=size(y0);nt=size(times);call init_real(sol,n,nt,times,y0)
    if(size(yprime0)/=n)then;sol%status=-101;sol%message='yprime0 size mismatch';return;end if
    if(.not.valid_times(times))then;sol%status=-100;sol%message='times must be strictly monotone';return;end if
    lrw=max(500,200+50*n+4*n*n);liw=max(100,100+10*n)
    allocate(rwork(lrw),iwork(liw),y(n),yp(n));rwork=0.0_dp;iwork=0;info=0;y=y0;yp=yprime0
    ! deSolve extends DDASPK INFO to 25 entries; INFO(21:23) stores nind.
    ! The ordinary ODE/DAE default is nind = [n, 0, 0].
    info(21)=n;info(22)=0;info(23)=0
    info(4)=1
    rt=1e-8_dp;at=1e-10_dp;if(present(rtol))rt=rtol;if(present(atol))at=atol
    idid=1;t=times(1);rpar=0.0_dp;ipar=0;current_residual=>residual;current_n=n
    do i=2,nt
      if(i==2)then;info(1)=0;else;info(1)=1;end if
      call ddaspk(daspk_res_bridge,n,t,y,yp,times(i),info,rt,at,idid,rwork,lrw,iwork,liw, &
          rpar,ipar,daspk_jac_dummy,daspk_psol_dummy)
      sol%y(:,i)=y
      if(idid<0)exit
    end do
    nullify(current_residual);current_n=0;sol%status=idid
    if(idid>0)then;sol%message='success';else;sol%message='DDASPK failed with IDID='//itoa(idid);end if
    if(liw>=13)then;sol%stats%n_steps=iwork(11);sol%stats%n_rhs=iwork(12);sol%stats%n_jac=iwork(13);end if
    if(lrw>=12)then;sol%stats%step_last=rwork(11);sol%stats%step_next=rwork(12);end if
  end function daspk

  function zvode(rhs,y0,times,rtol,atol,mf) result(sol)
    procedure(complex_ode_rhs)::rhs
    complex(dp),intent(in)::y0(:)
    real(dp),intent(in)::times(:)
    real(dp),intent(in),optional::rtol,atol
    integer,intent(in),optional::mf
    type(complex_ode_result)::sol
    integer::n,nt,i,itol,itask,istate,iopt,lzw,lrw,liw,mf0
    integer,allocatable::iwork(:)
    complex(dp),allocatable::zwork(:),y(:)
    real(dp),allocatable::rwork(:)
    real(dp)::rt(1),at(1),t,rpar(1)
    integer::ipar(1)
    interface
      subroutine desolve_zvode_entry(f,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,zwork,lzw, &
          rwork,lrw,iwork,liw,jac,mf,rpar,ipar)
        import dp
        external::f,jac
        integer::neq,itol,itask,istate,iopt,lzw,lrw,iwork(*),liw,mf,ipar(*)
        complex(dp)::y(*),zwork(*)
        real(dp)::t,tout,rtol(*),atol(*),rwork(*),rpar(*)
      end subroutine desolve_zvode_entry
    end interface
    n=size(y0);nt=size(times);allocate(sol%t(nt),sol%y(n,nt));sol%t=times;sol%y=(0.0_dp,0.0_dp);sol%y(:,1)=y0
    if(.not.valid_times(times))then;sol%status=-100;sol%message='times must be strictly monotone';return;end if
    itol=1;itask=1;istate=1;iopt=0;mf0=22;if(present(mf))mf0=mf
    lzw=max(100,100+20*n+n*n);lrw=max(100,100+20*n);liw=max(50,30+n)
    allocate(zwork(lzw),rwork(lrw),iwork(liw),y(n));zwork=(0.0_dp,0.0_dp);rwork=0.0_dp;iwork=0;y=y0
    rt=1e-8_dp;at=1e-10_dp;if(present(rtol))rt=rtol;if(present(atol))at=atol
    t=times(1);rpar=0.0_dp;ipar=0;current_complex_rhs=>rhs
    do i=2,nt
      call desolve_zvode_entry(zvode_rhs_bridge,n,y,t,times(i),itol,rt,at,itask,istate,iopt,zwork,lzw, &
          rwork,lrw,iwork,liw,zvode_jac_dummy,mf0,rpar,ipar)
      sol%y(:,i)=y;if(istate<0)exit
    end do
    nullify(current_complex_rhs);sol%status=istate
    if(istate==2)then;sol%message='success';else;sol%message='ZVODE status='//itoa(istate);end if
    if(liw>=14)then;sol%stats%n_steps=iwork(11);sol%stats%n_rhs=iwork(12);sol%stats%n_jac=iwork(13); &
        sol%stats%order_last=iwork(14);end if
    if(lrw>=12)then;sol%stats%step_last=rwork(11);sol%stats%step_next=rwork(12);end if
  end function zvode

  subroutine init_real(sol,n,nt,times,y0)
    type(ode_result),intent(out)::sol;integer,intent(in)::n,nt;real(dp),intent(in)::times(:),y0(:)
    allocate(sol%t(nt),sol%y(n,nt));sol%t=times;sol%y=0.0_dp;if(n>0.and.nt>0)sol%y(:,1)=y0
    sol%status=0;sol%message='not run'
  end subroutine init_real

  pure logical function valid_times(t) result(ok)
    real(dp),intent(in)::t(:);integer::i
    ok=size(t)>=1;if(size(t)<2)return
    if(t(2)>t(1))then
      do i=2,size(t);if(t(i)<=t(i-1))then;ok=.false.;return;end if;end do
    else if(t(2)<t(1))then
      do i=2,size(t);if(t(i)>=t(i-1))then;ok=.false.;return;end if;end do
    else;ok=.false.;end if
  end function valid_times

  pure function itoa(i) result(s)
    integer,intent(in)::i;character(len=:),allocatable::s;character(len=32)::b
    write(b,'(i0)')i;s=trim(b)
  end function itoa
end module desolve_stiff
