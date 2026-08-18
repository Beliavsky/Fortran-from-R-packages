! SPDX-License-Identifier: GPL-3.0-or-later
! Based on 'statnet' project software (statnet.org).
module degreenet_fit
  use degreenet_kinds, only : dp, huge_neg
  use degreenet_models, only : loglik_model
  implicit none
  private
  type, public :: fit_result
    real(dp), allocatable :: theta(:)
    real(dp), allocatable :: hessian(:,:), covariance(:,:), se(:)
    real(dp) :: loglik = huge_neg
    integer :: iterations = 0
    logical :: converged = .false.
  end type fit_result
  public :: fit_degree_model, nelder_mead_model

contains
  subroutine fit_degree_model(model,x,cutoff,cutabove,start,res,lower,upper,maxit,tol)
    integer,intent(in)::model,x(:),cutoff,cutabove
    real(dp),intent(in)::start(:)
    type(fit_result),intent(out)::res
    real(dp),intent(in),optional::lower(:),upper(:),tol
    integer,intent(in),optional::maxit
    integer::npar,i,j
    real(dp)::h,fp,fm,fpp,fpm,fmp,fmm,base
    npar=size(start);allocate(res%theta(npar),res%hessian(npar,npar), &
      res%covariance(npar,npar),res%se(npar))
    call nelder_mead_model(model,x,cutoff,cutabove,start,res%theta,res%loglik, &
      res%iterations,res%converged,lower,upper,maxit,tol)
    base=loglik_model(model,res%theta,x,cutoff,cutabove)
    res%hessian=0.0_dp
    do i=1,npar
      h=1e-4_dp*max(1.0_dp,abs(res%theta(i)))
      block
        real(dp)::t(npar)
        t=res%theta;t(i)=t(i)+h;fp=loglik_model(model,t,x,cutoff,cutabove)
        t=res%theta;t(i)=t(i)-h;fm=loglik_model(model,t,x,cutoff,cutabove)
      end block
      res%hessian(i,i)=(fp-2.0_dp*base+fm)/(h*h)
      do j=i+1,npar
        block
          real(dp)::t(npar),hj
          hj=1e-4_dp*max(1.0_dp,abs(res%theta(j)))
          t=res%theta;t(i)=t(i)+h;t(j)=t(j)+hj;fpp=loglik_model(model,t,x,cutoff,cutabove)
          t=res%theta;t(i)=t(i)+h;t(j)=t(j)-hj;fpm=loglik_model(model,t,x,cutoff,cutabove)
          t=res%theta;t(i)=t(i)-h;t(j)=t(j)+hj;fmp=loglik_model(model,t,x,cutoff,cutabove)
          t=res%theta;t(i)=t(i)-h;t(j)=t(j)-hj;fmm=loglik_model(model,t,x,cutoff,cutabove)
          res%hessian(i,j)=(fpp-fpm-fmp+fmm)/(4.0_dp*h*hj)
          res%hessian(j,i)=res%hessian(i,j)
        end block
      end do
    end do
    call inverse_matrix(-res%hessian,res%covariance,res%converged)
    do i=1,npar;res%se(i)=sqrt(max(0.0_dp,res%covariance(i,i)));end do
  end subroutine fit_degree_model

  subroutine nelder_mead_model(model,x,cutoff,cutabove,start,best,bestll,iters,converged, &
      lower,upper,maxit,tol)
    integer,intent(in)::model,x(:),cutoff,cutabove
    real(dp),intent(in)::start(:)
    real(dp),intent(out)::best(:),bestll
    integer,intent(out)::iters
    logical,intent(out)::converged
    real(dp),intent(in),optional::lower(:),upper(:),tol
    integer,intent(in),optional::maxit
    integer::n,m,j,ilo,ihi,inhi,mit
    real(dp)::ftol,alpha,gamma,rho,sigma,spread,fr,fe,fc
    real(dp),allocatable::simp(:,:),f(:),cent(:),xr(:),xe(:),xc(:)
    n=size(start);m=n+1;mit=3000;if(present(maxit))mit=maxit
    ftol=1e-8_dp;if(present(tol))ftol=tol
    alpha=1.0_dp;gamma=2.0_dp;rho=0.5_dp;sigma=0.5_dp
    allocate(simp(n,m),f(m),cent(n),xr(n),xe(n),xc(n));simp(:,1)=start
    do j=2,m
      simp(:,j)=start
      simp(j-1,j)=start(j-1)+0.08_dp*max(1.0_dp,abs(start(j-1)))
    end do
    do j=1,m;call clamp(simp(:,j));f(j)=-loglik_model(model,simp(:,j),x,cutoff,cutabove);end do
    converged=.false.
    do iters=1,mit
      ilo=minloc(f,dim=1);ihi=maxloc(f,dim=1);inhi=ilo
      do j=1,m;if(j/=ihi.and.(inhi==ihi.or.f(j)>f(inhi)))inhi=j;end do
      spread=maxval(abs(f-f(ilo)))/max(1.0_dp,abs(f(ilo)))
      if(spread<ftol)then;converged=.true.;exit;end if
      cent=0.0_dp;do j=1,m;if(j/=ihi)cent=cent+simp(:,j);end do;cent=cent/real(n,dp)
      xr=cent+alpha*(cent-simp(:,ihi));call clamp(xr);fr=-loglik_model(model,xr,x,cutoff,cutabove)
      if(fr<f(ilo))then
        xe=cent+gamma*(xr-cent);call clamp(xe);fe=-loglik_model(model,xe,x,cutoff,cutabove)
        if(fe<fr)then;simp(:,ihi)=xe;f(ihi)=fe;else;simp(:,ihi)=xr;f(ihi)=fr;end if
      else if(fr<f(inhi))then
        simp(:,ihi)=xr;f(ihi)=fr
      else
        if(fr<f(ihi))then;xc=cent+rho*(xr-cent);else;xc=cent+rho*(simp(:,ihi)-cent);end if
        call clamp(xc);fc=-loglik_model(model,xc,x,cutoff,cutabove)
        if(fc<min(fr,f(ihi)))then;simp(:,ihi)=xc;f(ihi)=fc
        else
          do j=1,m
            if(j/=ilo)then;simp(:,j)=simp(:,ilo)+sigma*(simp(:,j)-simp(:,ilo));call clamp(simp(:,j)); &
              f(j)=-loglik_model(model,simp(:,j),x,cutoff,cutabove);end if
          end do
        end if
      end if
    end do
    ilo=minloc(f,dim=1);best=simp(:,ilo);bestll=-f(ilo)
  contains
    subroutine clamp(v)
      real(dp),intent(inout)::v(:)
      if(present(lower))v=max(v,lower)
      if(present(upper))v=min(v,upper)
    end subroutine clamp
  end subroutine nelder_mead_model

  subroutine inverse_matrix(a,ainv,ok)
    real(dp),intent(in)::a(:,:);real(dp),intent(out)::ainv(:,:);logical,intent(out)::ok
    integer::n,i,j,k,p;real(dp)::mx,tmp
    real(dp),allocatable::aug(:,:),row(:)
    n=size(a,1);allocate(aug(n,2*n),row(2*n));aug=0.0_dp;aug(:,1:n)=a
    do i=1,n;aug(i,n+i)=1.0_dp;end do;ok=.true.
    do i=1,n
      p=i;mx=abs(aug(i,i));do k=i+1,n;if(abs(aug(k,i))>mx)then;mx=abs(aug(k,i));p=k;end if;end do
      if(mx<1e-12_dp)then;ok=.false.;ainv=0.0_dp;return;end if
      if(p/=i)then;row=aug(i,:);aug(i,:)=aug(p,:);aug(p,:)=row;end if
      tmp=aug(i,i);aug(i,:)=aug(i,:)/tmp
      do j=1,n;if(j/=i)then;tmp=aug(j,i);aug(j,:)=aug(j,:)-tmp*aug(i,:);end if;end do
    end do
    ainv=aug(:,n+1:2*n)
  end subroutine inverse_matrix
end module degreenet_fit
