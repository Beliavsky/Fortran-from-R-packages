! SPDX-License-Identifier: GPL-2.0-or-later
module tolerance_optimize
  use tolerance_kinds, only : dp
  implicit none
  private
  public :: golden_minimize, nelder_mead, numerical_hessian

  abstract interface
    function scalar_objective(x) result(v)
      import dp
      real(dp),intent(in)::x
      real(dp)::v
    end function scalar_objective
    function vector_objective(x) result(v)
      import dp
      real(dp),intent(in)::x(:)
      real(dp)::v
    end function vector_objective
  end interface

contains

  real(dp) function golden_minimize(f,a,b,tol) result(xmin)
    procedure(scalar_objective)::f
    real(dp),intent(in)::a,b
    real(dp),intent(in),optional::tol
    real(dp)::lo,hi,c,d,fc,fd,t,gr
    integer::iter
    t=1.0e-10_dp;if(present(tol))t=tol
    gr=(sqrt(5.0_dp)-1.0_dp)/2.0_dp;lo=a;hi=b
    c=hi-gr*(hi-lo);d=lo+gr*(hi-lo);fc=f(c);fd=f(d)
    do iter=1,400
      if(abs(hi-lo)<=t*max(1.0_dp,abs(0.5_dp*(lo+hi))))exit
      if(fc<fd)then
        hi=d;d=c;fd=fc;c=hi-gr*(hi-lo);fc=f(c)
      else
        lo=c;c=d;fc=fd;d=lo+gr*(hi-lo);fd=f(d)
      end if
    end do
    xmin=0.5_dp*(lo+hi)
  end function golden_minimize

  subroutine nelder_mead(f,x,step,tol,max_iter,fmin)
    procedure(vector_objective)::f
    real(dp),intent(inout)::x(:)
    real(dp),intent(in),optional::step,tol
    integer,intent(in),optional::max_iter
    real(dp),intent(out),optional::fmin
    real(dp),allocatable::simp(:,:),fv(:),cent(:),xr(:),xe(:),xc(:)
    real(dp)::st,tt,fr,fe,fc,tmpf,spread
    integer::n,i,j,k,mi
    real(dp),parameter::alpha=1.0_dp,gamma=2.0_dp,rho=0.5_dp,sigma=0.5_dp
    n=size(x);st=0.1_dp;if(present(step))st=step;tt=1.0e-9_dp;if(present(tol))tt=tol
    mi=1000;if(present(max_iter))mi=max_iter
    allocate(simp(n,n+1),fv(n+1),cent(n),xr(n),xe(n),xc(n))
    simp(:,1)=x
    do j=2,n+1
      simp(:,j)=x
      simp(j-1,j)=simp(j-1,j)+st*max(1.0_dp,abs(x(j-1)))
    end do
    do j=1,n+1;fv(j)=f(simp(:,j));end do
    do k=1,mi
      call sort_simplex(simp,fv)
      spread=maxval(abs(fv-fv(1)))
      if(spread<=tt*max(1.0_dp,abs(fv(1))))exit
      cent=sum(simp(:,1:n),dim=2)/real(n,dp)
      xr=cent+alpha*(cent-simp(:,n+1));fr=f(xr)
      if(fr<fv(1))then
        xe=cent+gamma*(xr-cent);fe=f(xe)
        if(fe<fr)then;simp(:,n+1)=xe;fv(n+1)=fe;else;simp(:,n+1)=xr;fv(n+1)=fr;end if
      else if(fr<fv(n))then
        simp(:,n+1)=xr;fv(n+1)=fr
      else
        if(fr<fv(n+1))then
          xc=cent+rho*(xr-cent)
        else
          xc=cent+rho*(simp(:,n+1)-cent)
        end if
        fc=f(xc)
        if(fc<min(fr,fv(n+1)))then
          simp(:,n+1)=xc;fv(n+1)=fc
        else
          do j=2,n+1
            simp(:,j)=simp(:,1)+sigma*(simp(:,j)-simp(:,1));fv(j)=f(simp(:,j))
          end do
        end if
      end if
    end do
    call sort_simplex(simp,fv);x=simp(:,1);if(present(fmin))fmin=fv(1)
  contains
    subroutine sort_simplex(s,vals)
      real(dp),intent(inout)::s(:,:),vals(:)
      real(dp)::tv
      real(dp),allocatable::col(:)
      integer::a,b
      allocate(col(size(s,1)))
      do a=2,size(vals)
        b=a
        do while(b>1)
          if(vals(b)>=vals(b-1)) exit
          tv=vals(b);vals(b)=vals(b-1);vals(b-1)=tv
          col=s(:,b);s(:,b)=s(:,b-1);s(:,b-1)=col;b=b-1
        end do
      end do
    end subroutine sort_simplex
  end subroutine nelder_mead

  subroutine numerical_hessian(f,x,hess,rel_step)
    procedure(vector_objective)::f
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::hess(:,:)
    real(dp),intent(in),optional::rel_step
    real(dp),allocatable::xp(:),xm(:),xpp(:),xpm(:),xmp(:),xmm(:)
    real(dp)::h,hi,hj,f0
    integer::n,i,j
    n=size(x);h=1.0e-5_dp;if(present(rel_step))h=rel_step
    allocate(xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n));hess=0.0_dp;f0=f(x)
    do i=1,n
      hi=h*max(1.0_dp,abs(x(i)));xp=x;xm=x;xp(i)=xp(i)+hi;xm(i)=xm(i)-hi
      hess(i,i)=(f(xp)-2.0_dp*f0+f(xm))/(hi*hi)
      do j=i+1,n
        hj=h*max(1.0_dp,abs(x(j)))
        xpp=x;xpm=x;xmp=x;xmm=x
        xpp(i)=xpp(i)+hi;xpp(j)=xpp(j)+hj
        xpm(i)=xpm(i)+hi;xpm(j)=xpm(j)-hj
        xmp(i)=xmp(i)-hi;xmp(j)=xmp(j)+hj
        xmm(i)=xmm(i)-hi;xmm(j)=xmm(j)-hj
        hess(i,j)=(f(xpp)-f(xpm)-f(xmp)+f(xmm))/(4.0_dp*hi*hj);hess(j,i)=hess(i,j)
      end do
    end do
  end subroutine numerical_hessian

end module tolerance_optimize
