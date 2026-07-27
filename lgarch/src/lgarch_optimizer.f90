! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2014-2025 Genaro Sucarrat (original R/C package)
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is distributed under the GNU General Public License version 2 only.
module lgarch_optimizer
  use lgarch_kinds, only : dp
  use lgarch_linalg, only : inverse_matrix
  implicit none
  private
  public :: minimize_nelder_mead, numerical_hessian, covariance_from_hessian

  abstract interface
    function objective_function(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function objective_function
  end interface
contains
  subroutine minimize_nelder_mead(fun,x,lower,upper,fmin,converged,iterations,max_iter,tol)
    procedure(objective_function) :: fun
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: lower(:),upper(:)
    real(dp), intent(out) :: fmin
    logical, intent(out) :: converged
    integer, intent(out) :: iterations
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tol
    integer :: n,mi,j,ilo,ihi,inhi
    real(dp) :: ftol,alpha,gamma,rho,sigma,spread,step
    real(dp), allocatable :: simplex(:,:),f(:),centroid(:),xr(:),xe(:),xc(:)
    n=size(x)
    if(size(lower)/=n .or. size(upper)/=n) error stop "minimize_nelder_mead: size mismatch"
    mi=5000; if(present(max_iter)) mi=max_iter
    ftol=1.0e-8_dp; if(present(tol)) ftol=tol
    alpha=1.0_dp; gamma=2.0_dp; rho=0.5_dp; sigma=0.5_dp
    allocate(simplex(n,n+1),f(n+1),centroid(n),xr(n),xe(n),xc(n))
    x=min(max(x,lower),upper); simplex(:,1)=x
    do j=1,n
      simplex(:,j+1)=x
      step=0.05_dp*max(1.0_dp,abs(x(j)))
      if(x(j)+step < upper(j)) then
        simplex(j,j+1)=x(j)+step
      else if(x(j)-step > lower(j)) then
        simplex(j,j+1)=x(j)-step
      else
        simplex(j,j+1)=0.5_dp*(lower(j)+upper(j))
      end if
    end do
    do j=1,n+1; f(j)=fun(simplex(:,j)); end do
    converged=.false.
    do iterations=1,mi
      call order_indices(f,ilo,ihi,inhi)
      spread=maxval(abs(f-f(ilo)))/(1.0_dp+abs(f(ilo)))
      if(spread<ftol .and. maxval(abs(simplex-spread_matrix(simplex(:,ilo),n+1)))<sqrt(ftol)*(1.0_dp+maxval(abs(simplex(:,ilo))))) then
        converged=.true.; exit
      end if
      centroid=0.0_dp
      do j=1,n+1; if(j/=ihi) centroid=centroid+simplex(:,j); end do
      centroid=centroid/real(n,dp)
      xr=clip(centroid+alpha*(centroid-simplex(:,ihi)),lower,upper)
      if(fun(xr)<f(ilo)) then
        xe=clip(centroid+gamma*(xr-centroid),lower,upper)
        if(fun(xe)<fun(xr)) then; simplex(:,ihi)=xe; f(ihi)=fun(xe); else; simplex(:,ihi)=xr; f(ihi)=fun(xr); end if
      else if(fun(xr)<f(inhi)) then
        simplex(:,ihi)=xr; f(ihi)=fun(xr)
      else
        if(fun(xr)<f(ihi)) then
          xc=clip(centroid+rho*(xr-centroid),lower,upper)
        else
          xc=clip(centroid-rho*(centroid-simplex(:,ihi)),lower,upper)
        end if
        if(fun(xc)<min(f(ihi),fun(xr))) then
          simplex(:,ihi)=xc; f(ihi)=fun(xc)
        else
          do j=1,n+1
            if(j/=ilo) then
              simplex(:,j)=clip(simplex(:,ilo)+sigma*(simplex(:,j)-simplex(:,ilo)),lower,upper)
              f(j)=fun(simplex(:,j))
            end if
          end do
        end if
      end if
    end do
    call order_indices(f,ilo,ihi,inhi); x=simplex(:,ilo); fmin=f(ilo)
  contains
    pure function clip(v,lo,hi) result(out)
      real(dp),intent(in)::v(:),lo(:),hi(:)
      real(dp)::out(size(v))
      out=min(max(v,lo),hi)
    end function clip
    pure function spread_matrix(v,m) result(a)
      real(dp),intent(in)::v(:); integer,intent(in)::m
      real(dp)::a(size(v),m); integer::k
      do k=1,m; a(:,k)=v; end do
    end function spread_matrix
    subroutine order_indices(vals,lo,hi,nhi)
      real(dp),intent(in)::vals(:); integer,intent(out)::lo,hi,nhi
      integer::k
      lo=minloc(vals,dim=1); hi=maxloc(vals,dim=1); nhi=merge(2,1,hi==1)
      do k=1,size(vals); if(k/=hi .and. vals(k)>vals(nhi)) nhi=k; end do
    end subroutine order_indices
  end subroutine minimize_nelder_mead

  subroutine numerical_hessian(fun,x,hess)
    procedure(objective_function) :: fun
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::hess(:,:)
    real(dp)::xp(size(x)),xm(size(x)),xpp(size(x)),xpm(size(x)),xmp(size(x)),xmm(size(x))
    real(dp)::hi,hj,f0
    integer::i,j,n
    n=size(x); if(any(shape(hess)/=[n,n])) error stop "numerical_hessian: size mismatch"
    f0=fun(x); hess=0.0_dp
    do i=1,n
      hi=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(i)))
      xp=x; xm=x; xp(i)=xp(i)+hi; xm(i)=xm(i)-hi
      hess(i,i)=(fun(xp)-2.0_dp*f0+fun(xm))/(hi*hi)
      do j=i+1,n
        hj=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(j)))
        xpp=x; xpm=x; xmp=x; xmm=x
        xpp(i)=xpp(i)+hi; xpp(j)=xpp(j)+hj
        xpm(i)=xpm(i)+hi; xpm(j)=xpm(j)-hj
        xmp(i)=xmp(i)-hi; xmp(j)=xmp(j)+hj
        xmm(i)=xmm(i)-hi; xmm(j)=xmm(j)-hj
        hess(i,j)=(fun(xpp)-fun(xpm)-fun(xmp)+fun(xmm))/(4.0_dp*hi*hj)
        hess(j,i)=hess(i,j)
      end do
    end do
  end subroutine numerical_hessian

  subroutine covariance_from_hessian(hess,cov,info)
    real(dp),intent(in)::hess(:,:)
    real(dp),intent(out)::cov(:,:)
    integer,intent(out)::info
    call inverse_matrix(hess,cov,info)
  end subroutine covariance_from_hessian
end module lgarch_optimizer
