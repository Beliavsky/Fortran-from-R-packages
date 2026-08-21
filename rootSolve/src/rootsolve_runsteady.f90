! SPDX-License-Identifier: GPL-2.0-or-later
module rootsolve_runsteady
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use rootsolve_kinds, only : dp
  use rootsolve_types, only : steady_rhs, runsteady_options, steady_result
  use desolve_callback_bridge, only : current_rhs, odepack_rhs_bridge, odepack_jac_dummy
  implicit none
  private
  public :: runsteady

  interface
    subroutine dlsode(f, neq, y, t, tout, itol, rtol, atol, itask, istate, &
        iopt, rwork, lrw, iwork, liw, jac, mf, rpar, ipar)
      import dp
      external :: f, jac
      integer :: neq(*), itol, itask, istate, iopt, lrw, iwork(*), liw, mf, ipar(*)
      real(dp) :: y(*), t, tout, rtol(*), atol(*), rwork(*), rpar(*)
    end subroutine dlsode
  end interface
contains

  function runsteady(func, y0, times, options) result(res)
    procedure(steady_rhs) :: func
    real(dp), intent(in) :: y0(:), times(2)
    type(runsteady_options), intent(in), optional :: options
    type(steady_result) :: res
    type(runsteady_options) :: opt
    real(dp), allocatable :: y(:),dy(:),rwork(:)
    integer, allocatable :: iwork(:)
    real(dp)::t,tout,rt(1),at(1),rpar(1),prec
    integer::n,neq(1),itol,itask,istate,iopt,lrw,liw,ipar(1),k

    if(present(options))opt=options
    n=size(y0)
    allocate(y(n),dy(n))
    y=y0
    neq(1)=n
    t=times(1)
    if(ieee_is_finite(times(2)))then
      tout=times(2)
    else
      if(times(2)>=times(1))then
        tout=huge(1.0_dp)/1024.0_dp
      else
        tout=-huge(1.0_dp)/1024.0_dp
      end if
    end if
    if(.not.(tout>t.or.tout<t))error stop 'runsteady: times must span a nonzero interval'

    itol=1
    itask=2
    istate=1
    iopt=1
    lrw=max(100,100+30*n+n*n)
    liw=max(50,30+n)
    allocate(rwork(lrw),iwork(liw),res%precision(opt%maxsteps))
    rwork=0.0_dp
    iwork=0
    res%precision=0.0_dp
    rt(1)=opt%rtol
    at(1)=opt%atol
    iwork(6)=opt%maxsteps
    rpar=0.0_dp
    ipar=0
    current_rhs=>func

    do k=1,opt%maxsteps
      call dlsode(odepack_rhs_bridge,neq,y,t,tout,itol,rt,at,itask,istate,iopt, &
          rwork,lrw,iwork,liw,odepack_jac_dummy,opt%mf,rpar,ipar)
      call func(t,y,dy)
      prec=sum(abs(dy))/real(max(1,n),dp)
      res%precision(k)=prec
      res%iterations=k
      if(prec<opt%stol)then
        res%steady=.true.
        exit
      end if
      if((tout>times(1).and.t>=tout).or.(tout<times(1).and.t<=tout))exit
      if(istate==-2)then
        rt=10.0_dp*rt
        at=10.0_dp*at
        istate=3
      else if(istate<0)then
        exit
      end if
      if(opt%positive)y=max(y,0.0_dp)
    end do
    nullify(current_rhs)

    allocate(res%y(n),res%f(n))
    res%y=y
    call func(t,y,res%f)
    res%time=t
    res%steps=iwork(11)
    if(res%iterations>0)res%estimated_precision=res%precision(res%iterations)
    if(res%iterations<size(res%precision))res%precision=res%precision(:res%iterations)
    if(res%steady)then
      res%status=1
    else if(istate<0)then
      res%status=istate
    else
      res%status=0
    end if
  end function runsteady
end module rootsolve_runsteady
