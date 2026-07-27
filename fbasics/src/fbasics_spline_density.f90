! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_spline_density
  use fbasics_kinds, only: dp, clamp
  use fbasics_rng, only: runif_lcg
  use fbasics_stats, only: sample_sd
  use fbasics_linalg, only: matrix_inverse
  implicit none
  private
  type, public :: spline_density_fit
    real(dp),allocatable::grid(:),density(:),cdf(:),coefficients(:),knots(:)
    real(dp)::lambda=0.0_dp,loglik=-huge(1.0_dp),penalized_loglik=-huge(1.0_dp)
    logical::converged=.false.
    integer::iterations=0
  end type spline_density_fit
  public :: fit_spline_density,dssd,pssd,qssd,rssd
contains
  subroutine fit_spline_density(x,fit,n_basis,lambda,grid_size,max_iter)
    real(dp),intent(in)::x(:)
    type(spline_density_fit),intent(out)::fit
    integer,intent(in),optional::n_basis,grid_size,max_iter
    real(dp),intent(in),optional::lambda
    integer::k,m,n,iter,imax,i,ls
    real(dp)::xmin,xmax,s,lam,stepnorm,oldobj,newobj
    real(dp),allocatable::bdata(:,:),bgrid(:,:),coef(:),grad(:),step(:),trial(:),p(:,:),info(:,:),invinfo(:,:),meanb(:),covb(:,:),dens(:),cdf(:)
    logical::accepted
    integer::info_inv
    n=size(x);k=12;if(present(n_basis))k=max(6,n_basis);m=501;if(present(grid_size))m=max(101,grid_size);imax=80;if(present(max_iter))imax=max_iter
    s=max(sample_sd(x),1.0e-6_dp);xmin=minval(x)-3.0_dp*s;xmax=maxval(x)+3.0_dp*s;if(xmax<=xmin)xmax=xmin+1.0_dp
    lam=0.1_dp*real(n,dp)**(-0.2_dp);if(present(lambda))lam=max(0.0_dp,lambda)
    call make_open_knots(xmin,xmax,k,3,fit%knots)
    allocate(fit%grid(m),bdata(n,k),bgrid(m,k),coef(k),grad(k),step(k),trial(k),p(k,k),info(k,k),meanb(k),covb(k,k),dens(m),cdf(m))
    do i=1,m;fit%grid(i)=xmin+(xmax-xmin)*real(i-1,dp)/real(m-1,dp);call bspline_basis(fit%grid(i),fit%knots,3,bgrid(i,:));end do
    do i=1,n;call bspline_basis(clamp(x(i),xmin,xmax),fit%knots,3,bdata(i,:));end do
    call second_difference_penalty(k,p);coef=0.0_dp;oldobj=penalized_objective(coef,bdata,bgrid,fit%grid,lam)
    fit%converged=.false.
    do iter=1,imax
      call density_moments(coef,bgrid,fit%grid,dens,meanb,covb)
      grad=sum(bdata,dim=1)-real(n,dp)*meanb-lam*matmul(p,coef)
      info=real(n,dp)*covb+lam*p
      do i=1,k;info(i,i)=info(i,i)+1.0e-8_dp;end do
      call matrix_inverse(info,invinfo,info_inv);if(info_inv/=0)exit;step=matmul(invinfo,grad);stepnorm=maxval(abs(step));accepted=.false.
      do ls=0,20
        trial=coef+step/(2.0_dp**ls);newobj=penalized_objective(trial,bdata,bgrid,fit%grid,lam)
        if(newobj>=oldobj)then;accepted=.true.;exit;end if
      end do
      if(.not.accepted)exit;coef=trial;oldobj=newobj
      if(stepnorm/(2.0_dp**ls)<1.0e-7_dp)then;fit%converged=.true.;exit;end if
    end do
    call density_moments(coef,bgrid,fit%grid,dens,meanb,covb);cdf(1)=0.0_dp
    do i=2,m;cdf(i)=cdf(i-1)+0.5_dp*(dens(i-1)+dens(i))*(fit%grid(i)-fit%grid(i-1));end do
    if(cdf(m)>0.0_dp)then;dens=dens/cdf(m);cdf=cdf/cdf(m);end if
    fit%density=dens;fit%cdf=cdf;fit%coefficients=coef;fit%lambda=lam;fit%iterations=iter
    fit%loglik=sum(matmul(bdata,coef))-real(n,dp)*log_partition(coef,bgrid,fit%grid)
    fit%penalized_loglik=fit%loglik-0.5_dp*lam*dot_product(coef,matmul(p,coef))
  end subroutine fit_spline_density

  subroutine make_open_knots(xmin,xmax,n_basis,degree,knots)
    real(dp),intent(in)::xmin,xmax;integer,intent(in)::n_basis,degree;real(dp),allocatable,intent(out)::knots(:)
    integer::nt,i,nint
    nt=n_basis+degree+1;allocate(knots(nt));knots(1:degree+1)=xmin;knots(nt-degree:nt)=xmax;nint=n_basis-degree-1
    do i=1,nint;knots(degree+1+i)=xmin+(xmax-xmin)*real(i,dp)/real(nint+1,dp);end do
  end subroutine make_open_knots

  subroutine bspline_basis(x,knots,degree,basis)
    real(dp),intent(in)::x,knots(:);integer,intent(in)::degree;real(dp),intent(out)::basis(:)
    real(dp),allocatable::prev(:),curr(:);real(dp)::left,right;integer::i,d,n
    n=size(basis);allocate(prev(n+degree),curr(n+degree));prev=0.0_dp
    do i=1,n+degree
      if((x>=knots(i).and.x<knots(i+1)).or.(x==knots(size(knots)).and.i==n))prev(i)=1.0_dp
    end do
    do d=1,degree
      curr=0.0_dp
      do i=1,n+degree-d
        left=0.0_dp;right=0.0_dp
        if(knots(i+d)>knots(i))left=(x-knots(i))/(knots(i+d)-knots(i))*prev(i)
        if(knots(i+d+1)>knots(i+1))right=(knots(i+d+1)-x)/(knots(i+d+1)-knots(i+1))*prev(i+1)
        curr(i)=left+right
      end do
      prev=curr
    end do
    basis=prev(1:n)
  end subroutine bspline_basis

  subroutine second_difference_penalty(k,p)
    integer,intent(in)::k;real(dp),intent(out)::p(k,k);real(dp),allocatable::d(:,:);integer::i
    allocate(d(k-2,k));d=0.0_dp;do i=1,k-2;d(i,i)=1.0_dp;d(i,i+1)=-2.0_dp;d(i,i+2)=1.0_dp;end do;p=matmul(transpose(d),d)
  end subroutine second_difference_penalty

  subroutine density_moments(coef,bgrid,grid,dens,meanb,covb)
    real(dp),intent(in)::coef(:),bgrid(:,:),grid(:);real(dp),intent(out)::dens(:),meanb(:),covb(:,:)
    real(dp),allocatable::eta(:),raw(:);real(dp)::z,w;integer::i,j,k,m
    m=size(grid);k=size(coef);allocate(eta(m),raw(m));eta=matmul(bgrid,coef);raw=exp(eta-maxval(eta));z=trapz(grid,raw);dens=raw/z;meanb=0.0_dp;covb=0.0_dp
    do i=2,m
      w=0.5_dp*(grid(i)-grid(i-1))
      meanb=meanb+w*(dens(i-1)*bgrid(i-1,:)+dens(i)*bgrid(i,:))
      do j=1,k;covb(j,:)=covb(j,:)+w*(dens(i-1)*bgrid(i-1,j)*bgrid(i-1,:)+dens(i)*bgrid(i,j)*bgrid(i,:));end do
    end do
    do j=1,k;covb(j,:)=covb(j,:)-meanb(j)*meanb;end do
  end subroutine density_moments

  real(dp) function log_partition(coef,bgrid,grid) result(v)
    real(dp),intent(in)::coef(:),bgrid(:,:),grid(:);real(dp),allocatable::eta(:);real(dp)::mx
    eta=matmul(bgrid,coef);mx=maxval(eta);v=mx+log(trapz(grid,exp(eta-mx)))
  end function log_partition
  real(dp) function penalized_objective(coef,bdata,bgrid,grid,lambda) result(v)
    real(dp),intent(in)::coef(:),bdata(:,:),bgrid(:,:),grid(:),lambda;real(dp),allocatable::p(:,:)
    allocate(p(size(coef),size(coef)));call second_difference_penalty(size(coef),p);v=sum(matmul(bdata,coef))-real(size(bdata,1),dp)*log_partition(coef,bgrid,grid)-0.5_dp*lambda*dot_product(coef,matmul(p,coef))
  end function penalized_objective
  pure real(dp) function trapz(x,y) result(v)
    real(dp),intent(in)::x(:),y(:);integer::i;v=0.0_dp;do i=2,size(x);v=v+0.5_dp*(y(i-1)+y(i))*(x(i)-x(i-1));end do
  end function trapz

  real(dp) function dssd(x,fit) result(v)
    real(dp),intent(in)::x;type(spline_density_fit),intent(in)::fit;v=linear_value(x,fit%grid,fit%density,0.0_dp,0.0_dp)
  end function dssd
  real(dp) function pssd(x,fit) result(v)
    real(dp),intent(in)::x;type(spline_density_fit),intent(in)::fit;v=linear_value(x,fit%grid,fit%cdf,0.0_dp,1.0_dp)
  end function pssd
  real(dp) function qssd(p,fit) result(v)
    real(dp),intent(in)::p;type(spline_density_fit),intent(in)::fit;integer::lo,hi,mid
    if(p<=0.0_dp)then;v=fit%grid(1);return;else if(p>=1.0_dp)then;v=fit%grid(size(fit%grid));return;end if
    lo=1;hi=size(fit%cdf);do while(hi-lo>1);mid=(lo+hi)/2;if(fit%cdf(mid)<p)then;lo=mid;else;hi=mid;end if;end do
    v=fit%grid(lo)+(fit%grid(hi)-fit%grid(lo))*(p-fit%cdf(lo))/max(fit%cdf(hi)-fit%cdf(lo),tiny(1.0_dp))
  end function qssd
  real(dp) function rssd(fit) result(v)
    type(spline_density_fit),intent(in)::fit;v=qssd(runif_lcg(),fit)
  end function rssd
  pure real(dp) function linear_value(x,g,y,left,right) result(v)
    real(dp),intent(in)::x,g(:),y(:),left,right;integer::lo,hi,mid
    if(x<=g(1))then;v=left;return;else if(x>=g(size(g)))then;v=right;return;end if;lo=1;hi=size(g);do while(hi-lo>1);mid=(lo+hi)/2;if(g(mid)<=x)then;lo=mid;else;hi=mid;end if;end do;v=y(lo)+(y(hi)-y(lo))*(x-g(lo))/(g(hi)-g(lo))
  end function linear_value
end module fbasics_spline_density
