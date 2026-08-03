! SPDX-License-Identifier: GPL-2.0-or-later
module jrvfinance_roots
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use jrvfinance_kinds, only: dp
  use jrvfinance_types, only: root_result, JRV_OK, JRV_INVALID_ARGUMENT, &
    JRV_NO_CONVERGENCE, JRV_NO_ROOT, JRV_NONFINITE
  implicit none
  private
  public :: root_callback, newton_raphson_root, bisection_root, irr_solve

  abstract interface
    subroutine root_callback(x, context, value, gradient)
      import dp
      real(dp), intent(in) :: x
      class(*), intent(in) :: context
      real(dp), intent(out) :: value, gradient
    end subroutine root_callback
  end interface
contains
  function newton_raphson_root(f, context, guess, lower, upper, max_iter, toler, convergence) result(res)
    procedure(root_callback) :: f
    class(*), intent(in) :: context
    real(dp), intent(in), optional :: guess, lower, upper, toler, convergence
    integer, intent(in), optional :: max_iter
    type(root_result) :: res
    real(dp) :: x, lo, hi, ftol, xtol, value, grad, step, trial, vtrial
    integer :: iter, nmax
    lo = -huge(1.0_dp); hi = huge(1.0_dp); x = 0.0_dp
    if (present(lower)) lo = lower
    if (present(upper)) hi = upper
    if (present(guess)) x = guess
    ftol = 1.0e-6_dp; xtol = 1.0e-8_dp; nmax = 100
    if (present(toler)) ftol = toler
    if (present(convergence)) xtol = convergence
    if (present(max_iter)) nmax = max_iter
    if (lo > x .or. hi < x .or. lo >= hi) then
      call fail(res, JRV_INVALID_ARGUMENT, 'guess must lie inside valid bounds')
      return
    end if
    do iter=0,nmax
      call f(x, context, value, grad)
      if (.not. ieee_is_finite(value) .or. .not. ieee_is_finite(grad)) then
        call fail(res, JRV_NONFINITE, 'callback returned a nonfinite value')
        return
      end if
      if (abs(value) <= ftol) then
        res%root=x; res%value=value; res%iterations=iter; res%status=JRV_OK
        return
      end if
      if (abs(grad) <= tiny(1.0_dp)) then
        trial = 0.5_dp*(lo+hi)
      else
        step = value/grad
        trial = x-step
        if (.not. ieee_is_finite(trial) .or. trial <= lo .or. trial >= hi) trial=0.5_dp*(lo+hi)
      end if
      call f(trial, context, vtrial, grad)
      if (.not. ieee_is_finite(vtrial)) then
        call fail(res, JRV_NONFINITE, 'callback returned a nonfinite value')
        return
      end if
      step = x - trial
      if (step > 0.0_dp) then
        hi = x
      else
        lo = x
      end if
      if (abs(trial-x) <= xtol*max(1.0_dp,abs(trial))) then
        x=trial
        call f(x, context, value, grad)
        if(abs(value)<=ftol) then
          res%root=x; res%value=value; res%iterations=iter+1; res%status=JRV_OK
        else
          call fail(res,JRV_NO_CONVERGENCE,'iteration stalled away from a root')
          res%root=x; res%value=value; res%iterations=iter+1
        end if
        return
      end if
      x=trial
    end do
    call f(x,context,value,grad)
    call fail(res,JRV_NO_CONVERGENCE,'maximum iterations exceeded')
    res%root=x; res%value=value; res%iterations=nmax
  end function newton_raphson_root

  function bisection_root(f, context, guess, lower, upper, nstep, toler) result(res)
    procedure(root_callback) :: f
    class(*), intent(in) :: context
    real(dp), intent(in) :: guess, lower, upper
    integer, intent(in), optional :: nstep
    real(dp), intent(in), optional :: toler
    type(root_result) :: res
    real(dp) :: left_ratio,right_ratio,low,up,l,r,fg,fl,fr,mid,fm,grad,ftol
    integer :: i,n
    logical :: bracketed
    n=100; if(present(nstep)) n=nstep
    ftol=1.0e-6_dp; if(present(toler)) ftol=toler
    if(lower<=0.0_dp .or. upper<=0.0_dp .or. guess<=0.0_dp .or. &
       lower>guess .or. upper<guess .or. n<1) then
      call fail(res,JRV_INVALID_ARGUMENT,'positive ordered bounds and guess are required')
      return
    end if
    call f(guess,context,fg,grad)
    if(abs(fg)<=ftol) then
      res%root=guess; res%value=fg; res%status=JRV_OK
      return
    end if
    left_ratio=(guess/lower)**(1.0_dp/real(n,dp))
    right_ratio=(upper/guess)**(1.0_dp/real(n,dp))
    low=guess; up=guess; l=0.0_dp; r=0.0_dp; bracketed=.false.
    do i=1,n
      low=low/left_ratio
      call f(low,context,fl,grad)
      if(opposite_sign(fl,fg)) then
        l=low; r=low*left_ratio; bracketed=.true.; exit
      end if
      up=up*right_ratio
      call f(up,context,fr,grad)
      if(opposite_sign(fr,fg)) then
        l=up/right_ratio; r=up; bracketed=.true.; exit
      end if
    end do
    if(.not.bracketed) then
      call fail(res,JRV_NO_ROOT,'failed to bracket a sign change')
      return
    end if
    call f(l,context,fl,grad)
    do i=1,200
      mid=0.5_dp*(l+r)
      call f(mid,context,fm,grad)
      if(abs(fm)<=ftol .or. abs(r-l)<=ftol*max(1.0_dp,abs(mid))) then
        res%root=mid; res%value=fm; res%iterations=i; res%status=JRV_OK
        return
      end if
      if(same_sign(fm,fl)) then
        l=mid; fl=fm
      else
        r=mid
      end if
    end do
    call fail(res,JRV_NO_CONVERGENCE,'bisection iteration limit exceeded')
    res%root=mid; res%value=fm; res%iterations=200
  end function bisection_root

  function irr_solve(f, context, interval, r_guess, toler, convergence, max_iter, method) result(res)
    procedure(root_callback) :: f
    class(*), intent(in) :: context
    real(dp), intent(in), optional :: interval(2), r_guess, toler, convergence
    integer, intent(in), optional :: max_iter
    character(len=*), intent(in), optional :: method
    type(root_result) :: res
    real(dp) :: lo,hi,guess,ftol,xtol,one,u_guess
    integer :: nmax
    character(len=16) :: how
    lo=-1.0_dp; hi=sqrt(real(huge(1),dp)); guess=0.0_dp
    if(present(interval)) then
      lo=max(interval(1),-1.0_dp); hi=interval(2); guess=0.5_dp*(lo+hi)
    end if
    if(present(r_guess)) guess=r_guess
    ftol=1.0e-6_dp; xtol=1.0e-8_dp; nmax=100; how='default'
    if(present(toler)) ftol=toler
    if(present(convergence)) xtol=convergence
    if(present(max_iter)) nmax=max_iter
    if(present(method)) how=lower(trim(method))
    if(lo>guess .or. hi<guess) then
      call fail(res,JRV_INVALID_ARGUMENT,'guess must lie inside interval')
      return
    end if
    if(trim(how)/='bisection') then
      res=newton_raphson_root(f,context,guess,lo,hi,nmax,ftol,xtol)
      if(res%status==JRV_OK .or. trim(how)=='newton') return
    end if
    one=1.01_dp
    u_guess=guess+one
    res=shifted_bisection(f,context,u_guess,lo+one,hi+one,ftol)
    if(res%status==JRV_OK) res%root=res%root-one
  end function irr_solve

  function shifted_bisection(f,context,guess,lower,upper,toler) result(res)
    procedure(root_callback) :: f
    class(*), intent(in) :: context
    real(dp), intent(in) :: guess,lower,upper,toler
    type(root_result) :: res
    real(dp) :: ratio_l,ratio_r,low,up,l,r,fg,fl,fr,mid,fm,g
    integer :: i
    logical :: bracketed
    call f(guess-1.01_dp,context,fg,g)
    ratio_l=(guess/lower)**0.01_dp
    ratio_r=(upper/guess)**0.01_dp
    low=guess; up=guess; l=0.0_dp; r=0.0_dp; bracketed=.false.
    do i=1,100
      low=low/ratio_l; call f(low-1.01_dp,context,fl,g)
      if(opposite_sign(fl,fg)) then; l=low; r=low*ratio_l; bracketed=.true.; exit; end if
      up=up*ratio_r; call f(up-1.01_dp,context,fr,g)
      if(opposite_sign(fr,fg)) then; l=up/ratio_r; r=up; bracketed=.true.; exit; end if
    end do
    if(.not.bracketed) then; call fail(res,JRV_NO_ROOT,'failed to bracket a sign change'); return; end if
    call f(l-1.01_dp,context,fl,g)
    do i=1,200
      mid=0.5_dp*(l+r); call f(mid-1.01_dp,context,fm,g)
      if(abs(fm)<=toler .or. abs(r-l)<=toler*max(1.0_dp,abs(mid))) then
        res%root=mid; res%value=fm; res%iterations=i; res%status=JRV_OK; return
      end if
      if(same_sign(fm,fl)) then; l=mid; fl=fm; else; r=mid; end if
    end do
    call fail(res,JRV_NO_CONVERGENCE,'bisection iteration limit exceeded')
  end function shifted_bisection

  pure logical function opposite_sign(a,b) result(value)
    real(dp), intent(in) :: a,b
    value = (a < 0.0_dp .and. b >= 0.0_dp) .or. (a >= 0.0_dp .and. b < 0.0_dp)
  end function opposite_sign

  pure logical function same_sign(a,b) result(value)
    real(dp), intent(in) :: a,b
    value = .not. opposite_sign(a,b)
  end function same_sign

  subroutine fail(res,status,message)
    type(root_result), intent(inout) :: res
    integer,intent(in)::status
    character(len=*),intent(in)::message
    res%status=status; res%message=message
  end subroutine fail

  pure function lower(text) result(out)
    character(len=*),intent(in)::text
    character(len=len(text))::out
    integer::i,c
    do i = 1, len(text)
      c = iachar(text(i:i))
      if (c >= 65 .and. c <= 90) then
        out(i:i) = achar(c+32)
      else
        out(i:i) = text(i:i)
      end if
    end do
  end function lower
end module jrvfinance_roots
