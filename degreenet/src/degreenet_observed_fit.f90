! SPDX-License-Identifier: GPL-3.0-or-later
! Based on 'statnet' project software (statnet.org).
module degreenet_observed_fit
  use degreenet_kinds, only : dp, huge_neg
  use degreenet_observation, only : grouped_loglik, rounded_loglik
  implicit none
  private
  type, public :: observed_fit_result
    real(dp), allocatable :: theta(:)
    real(dp) :: loglik = huge_neg
    integer :: iterations = 0
    logical :: converged = .false.
  end type observed_fit_result
  public :: fit_grouped_model, fit_rounded_model

contains
  subroutine fit_grouped_model(model,x,cutoff,start,res,lower,upper,maxit,tol)
    integer,intent(in)::model,x(:),cutoff
    real(dp),intent(in)::start(:)
    type(observed_fit_result),intent(out)::res
    real(dp),intent(in),optional::lower(:),upper(:),tol
    integer,intent(in),optional::maxit
    allocate(res%theta(size(start)))
    call nm_observed(1,model,x,cutoff,1000,start,res%theta,res%loglik,res%iterations, &
      res%converged,lower,upper,maxit,tol)
  end subroutine fit_grouped_model

  subroutine fit_rounded_model(model,x,cutoff,cutabove,start,res,lower,upper,maxit,tol)
    integer,intent(in)::model,x(:),cutoff,cutabove
    real(dp),intent(in)::start(:)
    type(observed_fit_result),intent(out)::res
    real(dp),intent(in),optional::lower(:),upper(:),tol
    integer,intent(in),optional::maxit
    allocate(res%theta(size(start)))
    call nm_observed(2,model,x,cutoff,cutabove,start,res%theta,res%loglik,res%iterations, &
      res%converged,lower,upper,maxit,tol)
  end subroutine fit_rounded_model

  subroutine nm_observed(mode,model,x,cutoff,cutabove,start,best,bestll,iters,converged, &
      lower,upper,maxit,tol)
    integer,intent(in)::mode,model,x(:),cutoff,cutabove
    real(dp),intent(in)::start(:)
    real(dp),intent(out)::best(:),bestll
    integer,intent(out)::iters
    logical,intent(out)::converged
    real(dp),intent(in),optional::lower(:),upper(:),tol
    integer,intent(in),optional::maxit
    integer::n,m,j,ilo,ihi,inhi,mit
    real(dp)::ftol,spread,fr,fe,fc
    real(dp),allocatable::simp(:,:),f(:),cent(:),xr(:),xe(:),xc(:)
    n=size(start);m=n+1;mit=3000;if(present(maxit))mit=maxit
    ftol=1e-8_dp;if(present(tol))ftol=tol
    allocate(simp(n,m),f(m),cent(n),xr(n),xe(n),xc(n));simp(:,1)=start
    do j=2,m
      simp(:,j)=start
      simp(j-1,j)=start(j-1)+0.08_dp*max(1.0_dp,abs(start(j-1)))
    end do
    do j=1,m;call clamp(simp(:,j));f(j)=-obj(simp(:,j));end do
    converged=.false.
    do iters=1,mit
      ilo=minloc(f,dim=1);ihi=maxloc(f,dim=1);inhi=ilo
      do j=1,m
        if(j/=ihi.and.(inhi==ihi.or.f(j)>f(inhi)))inhi=j
      end do
      spread=maxval(abs(f-f(ilo)))/max(1.0_dp,abs(f(ilo)))
      if(spread<ftol)then;converged=.true.;exit;end if
      cent=0.0_dp
      do j=1,m;if(j/=ihi)cent=cent+simp(:,j);end do
      cent=cent/real(n,dp)
      xr=cent+(cent-simp(:,ihi));call clamp(xr);fr=-obj(xr)
      if(fr<f(ilo))then
        xe=cent+2.0_dp*(xr-cent);call clamp(xe);fe=-obj(xe)
        if(fe<fr)then;simp(:,ihi)=xe;f(ihi)=fe;else;simp(:,ihi)=xr;f(ihi)=fr;end if
      else if(fr<f(inhi))then
        simp(:,ihi)=xr;f(ihi)=fr
      else
        if(fr<f(ihi))then;xc=cent+0.5_dp*(xr-cent);else;xc=cent+0.5_dp*(simp(:,ihi)-cent);end if
        call clamp(xc);fc=-obj(xc)
        if(fc<min(fr,f(ihi)))then;simp(:,ihi)=xc;f(ihi)=fc
        else
          do j=1,m
            if(j/=ilo)then
              simp(:,j)=simp(:,ilo)+0.5_dp*(simp(:,j)-simp(:,ilo));call clamp(simp(:,j));f(j)=-obj(simp(:,j))
            end if
          end do
        end if
      end if
    end do
    ilo=minloc(f,dim=1);best=simp(:,ilo);bestll=-f(ilo)
  contains
    real(dp) function obj(v) result(val)
      real(dp),intent(in)::v(:)
      if(mode==1)then;val=grouped_loglik(model,v,x,cutoff)
      else;val=rounded_loglik(model,v,x,cutoff,cutabove);end if
    end function obj
    subroutine clamp(v)
      real(dp),intent(inout)::v(:)
      if(present(lower))v=max(v,lower)
      if(present(upper))v=min(v,upper)
    end subroutine clamp
  end subroutine nm_observed
end module degreenet_observed_fit
