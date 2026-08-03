! SPDX-License-Identifier: GPL-3.0-only
module esback_optimizer
  use esback_kinds, only: dp, huge_penalty
  use esback_types, only: esback_ok, esback_optimization_failed
  implicit none
  private
  public :: objective_function, nelder_mead

  abstract interface
    function objective_function(x,context) result(f)
      import :: dp
      real(dp),intent(in)::x(:)
      class(*),intent(in)::context
      real(dp)::f
    end function objective_function
  end interface
contains
  subroutine nelder_mead(fun,context,x0,step,maxit,tol,xbest,fbest,status,iterations)
    procedure(objective_function)::fun
    class(*),intent(in)::context
    real(dp),intent(in)::x0(:)
    real(dp),intent(in),optional::step(:)
    integer,intent(in),optional::maxit
    real(dp),intent(in),optional::tol
    real(dp),allocatable,intent(out)::xbest(:)
    real(dp),intent(out)::fbest
    integer,intent(out)::status
    integer,intent(out),optional::iterations
    integer::n,mx,it,j,ilo,ihi,inhi
    real(dp)::tt,fr,fe,fc,spread,diam
    real(dp),allocatable::simp(:,:),fv(:),cent(:),xr(:),xe(:),xc(:),st(:)
    n=size(x0);mx=4000;if(present(maxit))mx=maxit;tt=1.0e-8_dp;if(present(tol))tt=tol
    allocate(simp(n,n+1),fv(n+1),cent(n),xr(n),xe(n),xc(n),st(n),xbest(n))
    if(present(step))then;st=max(abs(step),1.0e-4_dp);else;st=0.1_dp*max(abs(x0),1.0_dp);end if
    simp(:,1)=x0
    do j=1,n;simp(:,j+1)=x0;simp(j,j+1)=simp(j,j+1)+st(j);end do
    do j=1,n+1;fv(j)=fun(simp(:,j),context);if(.not.(fv(j)<huge_penalty))fv(j)=huge_penalty;end do
    do it=1,mx
      ilo=minloc(fv,dim=1);ihi=maxloc(fv,dim=1)
      inhi=ilo
      do j=1,n+1
        if(j/=ihi .and. (inhi==ilo .or. fv(j)>fv(inhi)))inhi=j
      end do
      spread=maxval(abs(fv-fv(ilo)))/max(1.0_dp,abs(fv(ilo)))
      diam=0.0_dp
      do j=1,n+1;diam=max(diam,maxval(abs(simp(:,j)-simp(:,ilo))));end do
      if(spread<tt .and. diam<sqrt(tt)*(1.0_dp+maxval(abs(simp(:,ilo)))))exit
      cent=0.0_dp
      do j=1,n+1;if(j/=ihi)cent=cent+simp(:,j);end do
      cent=cent/real(n,dp)
      xr=cent+(cent-simp(:,ihi));fr=fun(xr,context)
      if(fr<fv(ilo))then
        xe=cent+2.0_dp*(xr-cent);fe=fun(xe,context)
        if(fe<fr)then;simp(:,ihi)=xe;fv(ihi)=fe;else;simp(:,ihi)=xr;fv(ihi)=fr;end if
      else if(fr<fv(inhi))then
        simp(:,ihi)=xr;fv(ihi)=fr
      else
        if(fr<fv(ihi))then;xc=cent+0.5_dp*(xr-cent);else;xc=cent+0.5_dp*(simp(:,ihi)-cent);end if
        fc=fun(xc,context)
        if(fc<min(fr,fv(ihi)))then;simp(:,ihi)=xc;fv(ihi)=fc
        else
          do j=1,n+1
            if(j/=ilo)then;simp(:,j)=simp(:,ilo)+0.5_dp*(simp(:,j)-simp(:,ilo));fv(j)=fun(simp(:,j),context);end if
          end do
        end if
      end if
    end do
    ilo=minloc(fv,dim=1);xbest=simp(:,ilo);fbest=fv(ilo)
    if(it<=mx .and. fbest<huge_penalty)then;status=esback_ok;else;status=esback_optimization_failed;end if
    if(present(iterations))iterations=min(it,mx)
  end subroutine nelder_mead
end module esback_optimizer
