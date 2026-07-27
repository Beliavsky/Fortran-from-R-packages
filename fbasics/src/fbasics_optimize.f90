! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_optimize
  use fbasics_kinds, only: dp, clamp
  implicit none
  private
  public :: nelder_mead_bounded, numerical_hessian
  abstract interface
    function vector_objective(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function vector_objective
  end interface
contains
  subroutine nelder_mead_bounded(fun,start,lower,upper,best,fbest,converged,max_iter,tol)
    procedure(vector_objective)::fun
    real(dp),intent(in)::start(:),lower(:),upper(:)
    real(dp),allocatable,intent(out)::best(:)
    real(dp),intent(out)::fbest
    logical,intent(out)::converged
    integer,intent(in),optional::max_iter
    real(dp),intent(in),optional::tol
    integer::n,j,it,imax
    real(dp)::alpha,gamma,rho,sigma,eps,fr,fe,fc
    real(dp),allocatable::simp(:,:),fv(:),cent(:),xr(:),xe(:),xc(:),tmp(:)
    n=size(start);imax=2000;if(present(max_iter))imax=max_iter;eps=1e-8_dp;if(present(tol))eps=tol
    alpha=1;gamma=2;rho=0.5;sigma=0.5
    allocate(simp(n,n+1),fv(n+1),cent(n),xr(n),xe(n),xc(n),tmp(n),best(n))
    simp(:,1)=min(max(start,lower),upper)
    do j=1,n;simp(:,j+1)=simp(:,1);simp(j,j+1)=clamp(simp(j,j+1)+0.05_dp*max(1.0_dp,abs(simp(j,j+1))),lower(j),upper(j));if(simp(j,j+1)==simp(j,1))simp(j,j+1)=clamp(simp(j,j+1)-0.05_dp,lower(j),upper(j));end do
    do j=1,n+1;fv(j)=fun(simp(:,j));end do
    converged=.false.
    do it=1,imax
      call sort_simplex(simp,fv)
      if(maxval(abs(fv-fv(1)))<eps.and.maxval(abs(simp-spread(simp(:,1),2,n+1)))<sqrt(eps))then;converged=.true.;exit;end if
      cent=sum(simp(:,1:n),dim=2)/real(n,dp);xr=min(max(cent+alpha*(cent-simp(:,n+1)),lower),upper);fr=fun(xr)
      if(fr<fv(1))then
        xe=min(max(cent+gamma*(xr-cent),lower),upper);fe=fun(xe)
        if(fe<fr)then;simp(:,n+1)=xe;fv(n+1)=fe;else;simp(:,n+1)=xr;fv(n+1)=fr;end if
      else if(fr<fv(n))then;simp(:,n+1)=xr;fv(n+1)=fr
      else
        if(fr<fv(n+1))then;xc=min(max(cent+rho*(xr-cent),lower),upper);else;xc=min(max(cent-rho*(cent-simp(:,n+1)),lower),upper);end if
        fc=fun(xc)
        if(fc<min(fr,fv(n+1)))then;simp(:,n+1)=xc;fv(n+1)=fc
        else
          do j=2,n+1;simp(:,j)=min(max(simp(:,1)+sigma*(simp(:,j)-simp(:,1)),lower),upper);fv(j)=fun(simp(:,j));end do
        end if
      end if
    end do
    call sort_simplex(simp,fv);best=simp(:,1);fbest=fv(1)
  contains
    subroutine sort_simplex(s,f)
      real(dp),intent(inout)::s(:,:),f(:);integer::a,b;real(dp)::ft
      do a=2,size(f);ft=f(a);tmp=s(:,a);b=a-1;do while(b>=1);if(f(b)<=ft)exit;f(b+1)=f(b);s(:,b+1)=s(:,b);b=b-1;end do;f(b+1)=ft;s(:,b+1)=tmp;end do
    end subroutine
  end subroutine

  subroutine numerical_hessian(fun,x,hess,step)
    procedure(vector_objective)::fun
    real(dp),intent(in)::x(:);real(dp),allocatable,intent(out)::hess(:,:);real(dp),intent(in),optional::step
    integer::i,j,n;real(dp)::h,hi,hj,fpp,fpm,fmp,fmm,f0
    real(dp),allocatable::xp(:),xm(:)
    n=size(x);h=1e-4_dp;if(present(step))h=step;allocate(hess(n,n),xp(n),xm(n));f0=fun(x)
    do i=1,n;hi=h*max(1.0_dp,abs(x(i)));xp=x;xm=x;xp(i)=xp(i)+hi;xm(i)=xm(i)-hi;hess(i,i)=(fun(xp)-2*f0+fun(xm))/(hi*hi)
      do j=1,i-1;hj=h*max(1.0_dp,abs(x(j)));xp=x;xp(i)=xp(i)+hi;xp(j)=xp(j)+hj;fpp=fun(xp);xp=x;xp(i)=xp(i)+hi;xp(j)=xp(j)-hj;fpm=fun(xp);xp=x;xp(i)=xp(i)-hi;xp(j)=xp(j)+hj;fmp=fun(xp);xp=x;xp(i)=xp(i)-hi;xp(j)=xp(j)-hj;fmm=fun(xp);hess(i,j)=(fpp-fpm-fmp+fmm)/(4*hi*hj);hess(j,i)=hess(i,j);end do
    end do
  end subroutine
end module fbasics_optimize
