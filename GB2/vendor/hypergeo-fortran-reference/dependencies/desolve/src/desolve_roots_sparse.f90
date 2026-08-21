! SPDX-License-Identifier: GPL-2.0-or-later
module desolve_roots_sparse
  use desolve_kinds, only : dp
  use desolve_types, only : ode_rhs, ode_root, ode_result, lsodar_result
  use desolve_callback_bridge, only : current_rhs, current_root, odepack_rhs_bridge, &
      odepack_jac_dummy, odepack_root_bridge, lsodes_jac_dummy
  implicit none
  private
  public :: lsodes, lsodar
  interface
    subroutine dlsodes(f,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw, &
        iwork,liw,iwk,jac,mf,rpar,ipar)
      import dp
      external::f,jac
      integer::neq(*),itol,itask,istate,iopt,lrw,iwork(*),liw,iwk(*),mf,ipar(*)
      real(dp)::y(*),t,tout,rtol(*),atol(*),rwork(*),rpar(*)
    end subroutine dlsodes
    subroutine dlsodar(f,neq,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw, &
        iwork,liw,jac,jt,g,ng,jroot,rpar,ipar)
      import dp
      external::f,jac,g
      integer::neq(*),itol,itask,istate,iopt,lrw,iwork(*),liw,jt,ng,jroot(*),ipar(*)
      real(dp)::y(*),t,tout,rtol(*),atol(*),rwork(*),rpar(*)
    end subroutine dlsodar
  end interface
contains
  function lsodes(rhs,y0,times,rtol,atol,mf,max_steps) result(sol)
    procedure(ode_rhs)::rhs
    real(dp),intent(in)::y0(:),times(:)
    real(dp),intent(in),optional::rtol,atol
    integer,intent(in),optional::mf,max_steps
    type(ode_result)::sol
    integer::n,nt,neq(1),itol,itask,istate,iopt,lrw,liw,mf0,i
    integer,allocatable::iwork(:),iwk(:)
    real(dp),allocatable::rwork(:),y(:)
    real(dp)::rt(1),at(1),t,rpar(1);integer::ipar(1)
    n=size(y0);nt=size(times);call init_sol(sol,n,nt,times,y0)
    if(.not.valid_times(times))then;sol%status=-100;sol%message='times must be strictly monotone';return;end if
    neq(1)=n;itol=1;itask=1;istate=1;iopt=0;mf0=222;if(present(mf))mf0=mf
    lrw=max(500,100+4*n*n+50*n);liw=max(200,100+2*n*n+50*n)
    allocate(rwork(lrw),iwork(liw),iwk(2*lrw),y(n));rwork=0.0_dp;iwork=0;iwk=0;y=y0
    rt=1e-8_dp;at=1e-10_dp;if(present(rtol))rt=rtol;if(present(atol))at=atol
    if(present(max_steps))then;iopt=1;iwork(6)=max_steps;end if
    t=times(1);rpar=0.0_dp;ipar=0;current_rhs=>rhs
    do i=2,nt
      call dlsodes(odepack_rhs_bridge,neq,y,t,times(i),itol,rt,at,itask,istate,iopt,rwork,lrw, &
          iwork,liw,iwk,lsodes_jac_dummy,mf0,rpar,ipar)
      sol%y(:,i)=y;if(istate<0)exit
    end do
    nullify(current_rhs);sol%status=istate
    if(istate==2)then;sol%message='success';else;sol%message='DLSODES status='//itoa(istate);end if
    call stats(sol,iwork,rwork)
  end function lsodes

  function lsodar(rhs,root,y0,times,ng,rtol,atol,max_roots) result(out)
    procedure(ode_rhs)::rhs
    procedure(ode_root)::root
    real(dp),intent(in)::y0(:),times(:)
    integer,intent(in)::ng
    real(dp),intent(in),optional::rtol,atol
    integer,intent(in),optional::max_roots
    type(lsodar_result)::out
    integer::n,nt,neq(1),itol,itask,istate,iopt,lrw,liw,jt,i,j,mr,nr
    integer,allocatable::iwork(:),jroot(:),ri(:)
    real(dp),allocatable::rwork(:),y(:),rtm(:),rst(:,:)
    real(dp)::rt(1),at(1),t,rpar(1);integer::ipar(1)
    n=size(y0);nt=size(times);call init_sol(out%solution,n,nt,times,y0)
    if(ng<1)then;out%solution%status=-102;out%solution%message='ng must be positive';return;end if
    if(.not.valid_times(times))then;out%solution%status=-100;out%solution%message='times must be strictly monotone';return;end if
    mr=100;if(present(max_roots))mr=max_roots
    allocate(rtm(mr),ri(mr),rst(n,mr));rtm=0.0_dp;ri=0;rst=0.0_dp;nr=0
    neq(1)=n;itol=1;itask=1;istate=1;iopt=0;jt=2
    lrw=max(500,100+4*n*n+50*n+3*ng);liw=max(100,50+n+ng)
    allocate(rwork(lrw),iwork(liw),jroot(ng),y(n));rwork=0.0_dp;iwork=0;jroot=0;y=y0
    rt=1e-8_dp;at=1e-10_dp;if(present(rtol))rt=rtol;if(present(atol))at=atol
    t=times(1);rpar=0.0_dp;ipar=0;current_rhs=>rhs;current_root=>root
    do i=2,nt
      do
        call dlsodar(odepack_rhs_bridge,neq,y,t,times(i),itol,rt,at,itask,istate,iopt,rwork,lrw, &
            iwork,liw,odepack_jac_dummy,jt,odepack_root_bridge,ng,jroot,rpar,ipar)
        if(istate==3)then
          do j=1,ng
            if(jroot(j)/=0.and.nr<mr)then;nr=nr+1;rtm(nr)=t;ri(nr)=j;rst(:,nr)=y;end if
          end do
          istate=2
        else
          exit
        end if
      end do
      out%solution%y(:,i)=y;if(istate<0)exit
    end do
    nullify(current_rhs);nullify(current_root);out%solution%status=istate
    if(istate==2)then;out%solution%message='success';else;out%solution%message='DLSODAR status='//itoa(istate);end if
    out%nroots=nr;allocate(out%root_time(nr),out%root_index(nr),out%root_state(n,nr))
    if(nr>0)then;out%root_time=rtm(:nr);out%root_index=ri(:nr);out%root_state=rst(:,:nr);end if
    call stats(out%solution,iwork,rwork)
  end function lsodar

  subroutine init_sol(sol,n,nt,times,y0)
    type(ode_result),intent(out)::sol;integer,intent(in)::n,nt;real(dp),intent(in)::times(:),y0(:)
    allocate(sol%t(nt),sol%y(n,nt));sol%t=times;sol%y=0.0_dp;if(n>0.and.nt>0)sol%y(:,1)=y0
    sol%status=0;sol%message='not run'
  end subroutine init_sol
  subroutine stats(sol,iw,rw)
    type(ode_result),intent(inout)::sol;integer,intent(in)::iw(:);real(dp),intent(in)::rw(:)
    if(size(iw)>=20)then;sol%stats%n_steps=iw(11);sol%stats%n_rhs=iw(12);sol%stats%n_jac=iw(13); &
      sol%stats%order_last=iw(14);sol%stats%method_last=iw(19);end if
    if(size(rw)>=15)then;sol%stats%step_last=rw(11);sol%stats%step_next=rw(12);sol%stats%time_last_switch=rw(15);end if
  end subroutine stats
  pure logical function valid_times(t) result(ok)
    real(dp),intent(in)::t(:);integer::i
    ok=size(t)>=1;if(size(t)<2)return
    if(t(2)>t(1))then;do i=2,size(t);if(t(i)<=t(i-1))then;ok=.false.;return;end if;end do
    else if(t(2)<t(1))then;do i=2,size(t);if(t(i)>=t(i-1))then;ok=.false.;return;end if;end do
    else;ok=.false.;end if
  end function valid_times
  pure function itoa(i) result(s)
    integer,intent(in)::i;character(len=:),allocatable::s;character(len=32)::b;write(b,'(i0)')i;s=trim(b)
  end function itoa
end module desolve_roots_sparse
