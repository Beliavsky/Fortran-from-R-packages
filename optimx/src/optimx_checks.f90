! SPDX-License-Identifier: GPL-2.0-only
module optimx_checks
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use optimx_kinds, only: dp
  use optimx_types
  use optimx_linalg, only: vector_norm, cholesky_pd, symmetric_min_bound
  use optimx_eval
  implicit none
  private
  public :: fnchk, grchk, hesschk, gHgen, gHgenb, pd_check, kktchk
  public :: scalechk, bmchk, optchk

contains
  subroutine fnchk(problem,x,value,ok,status)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::value
    logical,intent(out)::ok
    integer,intent(out)::status
    integer::count
    count=0
    if(.not.valid_problem(problem,x))then
      value=huge(1.0_dp);ok=.false.;status=OPTIMX_INVALID_INPUT;return
    end if
    call evaluate_value(problem,x,value,count,status)
    ok=status==0 .and. ieee_is_finite(value)
  end subroutine fnchk

  subroutine grchk(problem,x,check,tol)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::x(:)
    type(derivative_check),intent(out)::check
    real(dp),intent(in),optional::tol
    real(dp)::threshold
    integer::fc,gc,status
    threshold=1.0e-5_dp;if(present(tol))threshold=tol
    allocate(check%analytic(size(x)),check%numeric(size(x)),check%error(size(x)))
    fc=0;gc=0
    if(.not.problem%has_gradient)then
      check%status=OPTIMX_INVALID_INPUT;check%ok=.false.;return
    end if
    call evaluate_gradient(problem,x,check%analytic,fc,gc,status)
    if(status/=0)then;check%status=status;return;end if
    call grcentral(problem,x,check%numeric,fc,status)
    if(status/=0)then;check%status=status;return;end if
    check%error=abs(check%analytic-check%numeric)/max(1.0_dp,abs(check%analytic),abs(check%numeric))
    check%max_error=maxval(check%error);check%ok=check%max_error<=threshold;check%status=0
  end subroutine grchk

  subroutine hesschk(problem,x,check,tol)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::x(:)
    type(hessian_check),intent(out)::check
    real(dp),intent(in),optional::tol
    real(dp)::threshold,f,g(size(x))
    integer::fc,gc,hc,status
    threshold=1.0e-4_dp;if(present(tol))threshold=tol
    allocate(check%analytic(size(x),size(x)),check%numeric(size(x),size(x)),check%error(size(x),size(x)))
    if(.not.problem%has_hessian)then;check%status=OPTIMX_INVALID_INPUT;return;end if
    fc=0;gc=0;hc=0;g=0.0_dp;check%analytic=0.0_dp
    call problem%objective(x,f,g,check%analytic,.true.,.true.,status)
    if(status/=0)then;check%status=status;return;end if
    block
      type(optimx_problem)::numeric_problem
      numeric_problem=problem
      numeric_problem%has_hessian=.false.
      call evaluate_hessian(numeric_problem,x,check%numeric,fc,gc,hc,status)
    end block
    if(status/=0)then;check%status=status;return;end if
    check%error=abs(check%analytic-check%numeric)/max(1.0_dp,abs(check%analytic),abs(check%numeric))
    check%max_error=maxval(check%error);check%ok=check%max_error<=threshold;check%status=0
  end subroutine hesschk

  subroutine gHgen(problem,x,g,h,status)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::g(:),h(:,:)
    integer,intent(out)::status
    integer::fc,gc,hc
    fc=0;gc=0;hc=0
    call evaluate_gradient(problem,x,g,fc,gc,status);if(status/=0)return
    call evaluate_hessian(problem,x,h,fc,gc,hc,status)
  end subroutine gHgen

  subroutine gHgenb(problem,x,g,h,status)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::g(:),h(:,:)
    integer,intent(out)::status
    integer::i
    call gHgen(problem,x,g,h,status);if(status/=0)return
    do i=1,size(x)
      if(problem%mask(i)==0)then;g(i)=0.0_dp;h(i,:)=0.0_dp;h(:,i)=0.0_dp;end if
    end do
  end subroutine gHgenb

  subroutine pd_check(a,positive,minimum_bound,tol)
    real(dp),intent(in)::a(:,:)
    logical,intent(out)::positive
    real(dp),intent(out),optional::minimum_bound
    real(dp),intent(in),optional::tol
    real(dp),allocatable::l(:,:)
    real(dp)::threshold,bound
    integer::status,n
    n=size(a,1);threshold=1.0e-7_dp;if(present(tol))threshold=tol
    bound=symmetric_min_bound(0.5_dp*(a+transpose(a)))
    if(present(minimum_bound))minimum_bound=bound
    if(size(a,2)/=n)then;positive=.false.;return;end if
    allocate(l(n,n));call cholesky_pd(0.5_dp*(a+transpose(a)),l,status,threshold)
    positive=status==0
  end subroutine pd_check

  subroutine kktchk(problem,x,result,control)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::x(:)
    type(kkt_result),intent(out)::result
    type(optimx_control),intent(in),optional::control
    type(optimx_control)::ctrl
    real(dp)::g(size(x)),pg(size(x)),h(size(x),size(x))
    logical::pd
    integer::fc,gc,hc,status
    ctrl=optimx_control();if(present(control))ctrl=control;fc=0;gc=0;hc=0
    call evaluate_gradient(problem,x,g,fc,gc,status);if(status/=0)then;result%status=status;return;end if
    call projected_gradient(problem,x,g,pg);result%projected_gradient_norm=vector_norm(pg)
    result%kkt1=result%projected_gradient_norm<=ctrl%gradtol
    call evaluate_hessian(problem,x,h,fc,gc,hc,status);if(status/=0)then;result%status=status;return;end if
    call pd_check(h,pd,result%minimum_eigen_bound,100.0_dp*ctrl%gradtol);result%kkt2=pd;result%status=0
  end subroutine kktchk

  subroutine scalechk(par,lower,upper,result)
    real(dp),intent(in)::par(:),lower(:),upper(:)
    type(scale_result),intent(out)::result
    real(dp)::pmin,pmax,bmin,bmax
    logical::pany,bany
    integer::i
    pmin=huge(1.0_dp);pmax=0.0_dp;pany=.false.
    bmin=huge(1.0_dp);bmax=0.0_dp;bany=.false.
    do i=1,size(par)
      if(abs(par(i))>0.0_dp .and. ieee_is_finite(par(i)))then
        pmin=min(pmin,abs(par(i)));pmax=max(pmax,abs(par(i)));pany=.true.
      end if
      if(ieee_is_finite(lower(i)) .and. ieee_is_finite(upper(i)))then
        if(upper(i)>lower(i))then;bmin=min(bmin,upper(i)-lower(i));bmax=max(bmax,upper(i)-lower(i));bany=.true.;end if
      end if
    end do
    if(pany)result%parameter_ratio=pmax/pmin
    if(bany)result%bound_ratio=bmax/bmin
    result%well_scaled=max(result%parameter_ratio,result%bound_ratio)<=1.0e6_dp
  end subroutine scalechk

  subroutine bmchk(problem,x,result,shift_to_bound)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::x(:)
    type(bounds_result),intent(out)::result
    logical,intent(in),optional::shift_to_bound
    logical::shift
    integer::i
    shift=.true.;if(present(shift_to_bound))shift=shift_to_bound
    allocate(result%par(size(x)),result%mask(size(x)));result%par=x;result%mask=problem%mask
    result%feasible=.true.
    do i=1,size(x)
      if(problem%lower(i)>problem%upper(i))then;result%feasible=.false.;cycle;end if
      if(result%mask(i)==0)cycle
      if(x(i)<problem%lower(i))then
        result%feasible=.false.;if(shift)result%par(i)=problem%lower(i)
      else if(x(i)>problem%upper(i))then
        result%feasible=.false.;if(shift)result%par(i)=problem%upper(i)
      end if
      if(abs(problem%upper(i)-problem%lower(i))<=epsilon(1.0_dp)*max(1.0_dp,abs(problem%lower(i))))then
        result%par(i)=problem%lower(i);result%mask(i)=0
      end if
    end do
    result%status=0
  end subroutine bmchk

  subroutine optchk(problem,x,result,control)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::x(:)
    type(kkt_result),intent(out)::result
    type(optimx_control),intent(in),optional::control
    call kktchk(problem,x,result,control)
  end subroutine optchk
end module optimx_checks
