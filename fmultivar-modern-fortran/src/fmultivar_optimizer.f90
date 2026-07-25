! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fMultivar authors and translation contributors.
! This file is part of fmultivar-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
module fmultivar_optimizer
  use fmultivar_kinds, only : dp
  implicit none
  private
  public :: optimizer_result, nelder_mead
  type :: optimizer_result
    real(dp), allocatable :: x(:)
    real(dp) :: value = huge(1.0_dp)
    integer :: iterations = 0
    integer :: evaluations = 0
    logical :: converged = .false.
  end type optimizer_result
  abstract interface
    function objective_function(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function objective_function
  end interface
contains
  function nelder_mead(fun,x0,step,tol,max_iter,lower,upper) result(res)
    procedure(objective_function)::fun
    real(dp),intent(in)::x0(:)
    real(dp),intent(in),optional::step(:),tol,lower(:),upper(:)
    integer,intent(in),optional::max_iter
    type(optimizer_result)::res
    real(dp),allocatable::simplex(:,:),fval(:),centroid(:),xr(:),xe(:),xc(:),st(:)
    real(dp)::alpha,gamma,rho,sigma,toluse,fr,fe,fc,spread
    integer::n,itmax,i,j
    n=size(x0);itmax=3000;if(present(max_iter))itmax=max_iter
    toluse=1.0e-7_dp;if(present(tol))toluse=tol
    alpha=1.0_dp;gamma=2.0_dp;rho=0.5_dp;sigma=0.5_dp
    allocate(simplex(n,n+1),fval(n+1),centroid(n),xr(n),xe(n),xc(n),st(n))
    st=0.05_dp*max(1.0_dp,abs(x0));if(present(step))st=step
    simplex(:,1)=x0;call apply_bounds(simplex(:,1),lower,upper)
    do j=2,n+1
      simplex(:,j)=simplex(:,1);simplex(j-1,j)=simplex(j-1,j)+st(j-1)
      call apply_bounds(simplex(:,j),lower,upper)
    end do
    do j=1,n+1;fval(j)=safe_eval(fun,simplex(:,j));end do
    res%evaluations=n+1
    do i=1,itmax
      call sort_simplex(simplex,fval)
      spread=maxval(abs(fval-fval(1)))
      if(spread<=toluse*(1.0_dp+abs(fval(1))) .and. &
         maxval(abs(simplex-spread_col(simplex(:,1),n+1)))<=sqrt(toluse)*(1.0_dp+maxval(abs(simplex(:,1)))))then
        res%converged=.true.;exit
      end if
      centroid=sum(simplex(:,1:n),dim=2)/real(n,dp)
      xr=centroid+alpha*(centroid-simplex(:,n+1));call apply_bounds(xr,lower,upper)
      fr=safe_eval(fun,xr);res%evaluations=res%evaluations+1
      if(fr<fval(1))then
        xe=centroid+gamma*(xr-centroid);call apply_bounds(xe,lower,upper)
        fe=safe_eval(fun,xe);res%evaluations=res%evaluations+1
        if(fe<fr)then;simplex(:,n+1)=xe;fval(n+1)=fe
        else;simplex(:,n+1)=xr;fval(n+1)=fr;end if
      else if(fr<fval(n))then
        simplex(:,n+1)=xr;fval(n+1)=fr
      else
        if(fr<fval(n+1))then;xc=centroid+rho*(xr-centroid)
        else;xc=centroid+rho*(simplex(:,n+1)-centroid);end if
        call apply_bounds(xc,lower,upper);fc=safe_eval(fun,xc);res%evaluations=res%evaluations+1
        if(fc<min(fr,fval(n+1)))then
          simplex(:,n+1)=xc;fval(n+1)=fc
        else
          do j=2,n+1
            simplex(:,j)=simplex(:,1)+sigma*(simplex(:,j)-simplex(:,1))
            call apply_bounds(simplex(:,j),lower,upper)
            fval(j)=safe_eval(fun,simplex(:,j))
          end do
          res%evaluations=res%evaluations+n
        end if
      end if
      res%iterations=i
    end do
    call sort_simplex(simplex,fval)
    allocate(res%x(n));res%x=simplex(:,1);res%value=fval(1)
    if(.not.res%converged)res%iterations=itmax
  contains
    function spread_col(v,m) result(a)
      real(dp),intent(in)::v(:)
      integer,intent(in)::m
      real(dp)::a(size(v),m)
      integer::k
      do k=1,m;a(:,k)=v;end do
    end function spread_col
  end function nelder_mead

  function safe_eval(fun,x) result(f)
    procedure(objective_function)::fun
    real(dp),intent(in)::x(:)
    real(dp)::f
    f=fun(x)
    if(.not.(f<huge(1.0_dp)))f=huge(1.0_dp)/100.0_dp
  end function safe_eval

  subroutine apply_bounds(x,lower,upper)
    real(dp),intent(inout)::x(:)
    real(dp),intent(in),optional::lower(:),upper(:)
    if(present(lower))x=max(x,lower)
    if(present(upper))x=min(x,upper)
  end subroutine apply_bounds

  subroutine sort_simplex(x,f)
    real(dp),intent(inout)::x(:,:),f(:)
    integer::i,j,k
    real(dp)::tf
    real(dp),allocatable::tx(:)
    allocate(tx(size(x,1)))
    do i=1,size(f)-1
      k=i
      do j=i+1,size(f);if(f(j)<f(k))k=j;end do
      if(k/=i)then
        tf=f(i);f(i)=f(k);f(k)=tf;tx=x(:,i);x(:,i)=x(:,k);x(:,k)=tx
      end if
    end do
  end subroutine sort_simplex
end module fmultivar_optimizer
