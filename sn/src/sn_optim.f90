! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
module sn_optim
  use sn_kinds, only : dp
  use sn_status, only : sn_ok, sn_invalid_argument, sn_no_convergence
  use sn_linalg, only : inverse_general
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  implicit none
  private

  abstract interface
    function objective_function(x) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: value
    end function objective_function
  end interface

  public :: nelder_mead, numerical_gradient, numerical_hessian

contains

  subroutine nelder_mead(f,x0,xbest,fbest,info,iterations,max_iter,tol,step)
    procedure(objective_function) :: f
    real(dp), intent(in) :: x0(:)
    real(dp), allocatable, intent(out) :: xbest(:)
    real(dp), intent(out) :: fbest
    integer, intent(out) :: info
    integer, intent(out), optional :: iterations
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tol,step
    real(dp), allocatable :: simplex(:,:),fv(:),centroid(:),xr(:),xe(:),xc(:),tmpv(:)
    real(dp) :: alpha,gamma,rho,sigma,eps,stp,fr,fe,fc,max_spread
    integer :: n,maxit,iter,j,ilo,ihi,inhi

    n=size(x0)
    if (n<1) then
      allocate(xbest(0))
      fbest=huge(1.0_dp)
      info=sn_invalid_argument
      return
    end if
    maxit=4000
    if (present(max_iter)) maxit=max_iter
    eps=1.0e-8_dp
    if (present(tol)) eps=tol
    stp=0.1_dp
    if (present(step)) stp=step
    alpha=1.0_dp; gamma=2.0_dp; rho=0.5_dp; sigma=0.5_dp
    allocate(simplex(n,n+1),fv(n+1),centroid(n),xr(n),xe(n),xc(n),tmpv(n))
    simplex(:,1)=x0
    do j=1,n
      simplex(:,j+1)=x0
      simplex(j,j+1)=simplex(j,j+1)+stp*max(1.0_dp,abs(x0(j)))
    end do
    do j=1,n+1
      fv(j)=safe_eval(f,simplex(:,j))
    end do

    do iter=1,maxit
      call order_indices(fv,ilo,ihi,inhi)
      max_spread=0.0_dp
      do j=1,n+1
        max_spread=max(max_spread,maxval(abs(simplex(:,j)-simplex(:,ilo))))
      end do
      if (maxval(abs(fv-fv(ilo))) <= eps*max(1.0_dp,abs(fv(ilo))) .and. &
          max_spread <= sqrt(eps)*max(1.0_dp,maxval(abs(simplex(:,ilo))))) exit
      centroid=(sum(simplex,dim=2)-simplex(:,ihi))/real(n,dp)
      xr=centroid+alpha*(centroid-simplex(:,ihi))
      fr=safe_eval(f,xr)
      if (fr < fv(ilo)) then
        xe=centroid+gamma*(xr-centroid)
        fe=safe_eval(f,xe)
        if (fe < fr) then
          simplex(:,ihi)=xe; fv(ihi)=fe
        else
          simplex(:,ihi)=xr; fv(ihi)=fr
        end if
      else if (fr < fv(inhi)) then
        simplex(:,ihi)=xr; fv(ihi)=fr
      else
        if (fr < fv(ihi)) then
          xc=centroid+rho*(xr-centroid)
        else
          xc=centroid-rho*(centroid-simplex(:,ihi))
        end if
        fc=safe_eval(f,xc)
        if (fc < min(fr,fv(ihi))) then
          simplex(:,ihi)=xc; fv(ihi)=fc
        else
          tmpv=simplex(:,ilo)
          do j=1,n+1
            if (j==ilo) cycle
            simplex(:,j)=tmpv+sigma*(simplex(:,j)-tmpv)
            fv(j)=safe_eval(f,simplex(:,j))
          end do
        end if
      end if
    end do
    call order_indices(fv,ilo,ihi,inhi)
    allocate(xbest(n))
    xbest=simplex(:,ilo)
    fbest=fv(ilo)
    if (iter>maxit) then
      info=sn_no_convergence
    else
      info=sn_ok
    end if
    if (present(iterations)) iterations=min(iter,maxit)
  contains
    real(dp) function safe_eval(fun,x) result(v)
      procedure(objective_function) :: fun
      real(dp), intent(in) :: x(:)
      v=fun(x)
      if (ieee_is_nan(v) .or. abs(v)>=huge(1.0_dp)) v=0.25_dp*huge(1.0_dp)
    end function safe_eval

    subroutine order_indices(vals,lo,hi,nhi)
      real(dp), intent(in) :: vals(:)
      integer, intent(out) :: lo,hi,nhi
      integer :: k
      lo=1; hi=1
      do k=2,size(vals)
        if (vals(k)<vals(lo)) lo=k
        if (vals(k)>vals(hi)) hi=k
      end do
      nhi=merge(2,1,hi==1)
      do k=1,size(vals)
        if (k==hi) cycle
        if (vals(k)>vals(nhi)) nhi=k
      end do
    end subroutine order_indices
  end subroutine nelder_mead

  subroutine numerical_gradient(f,x,g,step)
    procedure(objective_function) :: f
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: g(:)
    real(dp), intent(in), optional :: step
    real(dp), allocatable :: xp(:),xm(:)
    real(dp) :: h
    integer :: n,i
    n=size(x)
    allocate(g(n),xp(n),xm(n))
    do i=1,n
      h=sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(x(i)))
      if (present(step)) h=step*max(1.0_dp,abs(x(i)))
      xp=x; xm=x
      xp(i)=xp(i)+h; xm(i)=xm(i)-h
      g(i)=(f(xp)-f(xm))/(2.0_dp*h)
    end do
  end subroutine numerical_gradient

  subroutine numerical_hessian(f,x,hessian,covariance,info,step)
    procedure(objective_function) :: f
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: hessian(:,:),covariance(:,:)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: step
    real(dp), allocatable :: xpp(:),xpm(:),xmp(:),xmm(:)
    real(dp) :: hi,hj,f0
    integer :: n,i,j
    n=size(x)
    allocate(hessian(n,n),xpp(n),xpm(n),xmp(n),xmm(n))
    f0=f(x)
    do i=1,n
      hi=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(i)))
      if (present(step)) hi=step*max(1.0_dp,abs(x(i)))
      xpp=x; xmm=x
      xpp(i)=xpp(i)+hi; xmm(i)=xmm(i)-hi
      hessian(i,i)=(f(xpp)-2.0_dp*f0+f(xmm))/(hi*hi)
      do j=i+1,n
        hj=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(j)))
        if (present(step)) hj=step*max(1.0_dp,abs(x(j)))
        xpp=x; xpm=x; xmp=x; xmm=x
        xpp(i)=xpp(i)+hi; xpp(j)=xpp(j)+hj
        xpm(i)=xpm(i)+hi; xpm(j)=xpm(j)-hj
        xmp(i)=xmp(i)-hi; xmp(j)=xmp(j)+hj
        xmm(i)=xmm(i)-hi; xmm(j)=xmm(j)-hj
        hessian(i,j)=(f(xpp)-f(xpm)-f(xmp)+f(xmm))/(4.0_dp*hi*hj)
        hessian(j,i)=hessian(i,j)
      end do
    end do
    call inverse_general(hessian,covariance,info)
  end subroutine numerical_hessian

end module sn_optim
