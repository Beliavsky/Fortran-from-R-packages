module discretedists_estimators
   use discretedists_kinds, only : dp
   use discretedists_numerics, only : nelder_mead, objective_fun, logistic, logit, mean_sample
   use discretedists_distributions
   implicit none
   private

   integer, parameter :: FIT_HYPERPO=1,FIT_HYPERPO2=2,FIT_DBH=3,FIT_DLD=4,FIT_DIKUM=5
   integer, parameter :: FIT_POISXL=6,FIT_GGEO=7,FIT_DGEII=8,FIT_DPERKS=9,FIT_DSPA=10,FIT_DMOLBE=11

   type :: fit_context_t
      integer :: family=0
      real(dp), allocatable :: y(:)
   end type fit_context_t

   public :: loglik_hyperpo,loglik_hyperpo2,loglik_dbh,loglik_dld,loglik_dikum
   public :: loglik_poisxl,loglik_ggeo,loglik_dgeii,loglik_dperks,loglik_dspa,loglik_dmolbe
   public :: estim_mu_sigma_hyperpo,estim_mu_sigma_hyperpo2,estim_mu_dbh,estim_mu_dld
   public :: estim_mu_sigma_dikum,estim_mu_poisxl,estim_mu_sigma_ggeo,estim_mu_sigma_dgeii
   public :: estim_mu_sigma_dperks,estim_mu_sigma_dspa,estim_mu_sigma_dmolbe,estim_mu_sigma_compo

contains

   function loglik_hyperpo(logparam,x) result(ll)
      real(dp),intent(in)::logparam(2),x(:)
      real(dp)::ll
      ll=sum(dhyperpo(x,exp(logparam(1)),exp(logparam(2)),.true.))
   end function loglik_hyperpo

   function loglik_hyperpo2(logparam,x) result(ll)
      real(dp),intent(in)::logparam(2),x(:)
      real(dp)::ll
      ll=sum(dhyperpo2(x,exp(logparam(1)),exp(logparam(2)),.true.))
   end function loglik_hyperpo2

   function loglik_dbh(param,x) result(nll)
      real(dp),intent(in)::param,x(:)
      real(dp)::nll
      nll=-sum(ddbh(x,param,.true.))
   end function loglik_dbh

   function loglik_dld(param,x) result(nll)
      real(dp),intent(in)::param,x(:)
      real(dp)::nll
      nll=-sum(ddld(x,param,.true.))
   end function loglik_dld

   function loglik_dikum(param,x) result(ll)
      real(dp),intent(in)::param(2),x(:)
      real(dp)::ll
      ll=sum(ddikum(x,exp(param(1)),exp(param(2)),.true.))
   end function loglik_dikum

   function loglik_poisxl(param,x) result(ll)
      real(dp),intent(in)::param,x(:)
      real(dp)::ll
      ll=sum(dpoisxl(x,exp(param),.true.))
   end function loglik_poisxl

   function loglik_ggeo(param,x) result(ll)
      real(dp),intent(in)::param(2),x(:)
      real(dp)::ll
      ll=sum(dggeo(x,logistic(param(1)),exp(param(2)),.true.))
   end function loglik_ggeo

   function loglik_dgeii(param,x) result(ll)
      real(dp),intent(in)::param(2),x(:)
      real(dp)::ll
      ll=sum(ddgeii(x,logistic(param(1)),exp(param(2)),.true.))
   end function loglik_dgeii

   function loglik_dperks(logparam,x) result(ll)
      real(dp),intent(in)::logparam(2),x(:)
      real(dp)::ll
      ll=sum(ddperks(x,exp(logparam(1)),exp(logparam(2)),.true.))
   end function loglik_dperks

   function loglik_dspa(param,x) result(ll)
      real(dp),intent(in)::param(2),x(:)
      real(dp)::ll
      ll=sum(ddspa(x,exp(param(1)),logistic(param(2)),.true.))
   end function loglik_dspa

   function loglik_dmolbe(logparam,x) result(ll)
      real(dp),intent(in)::logparam(2),x(:)
      real(dp)::ll
      ll=sum(ddmolbe(x,exp(logparam(1)),exp(logparam(2)),.true.))
   end function loglik_dmolbe

   function fit_objective(theta,context) result(value)
      real(dp),intent(in)::theta(:)
      class(*),intent(in)::context
      real(dp)::value
      select type(ctx=>context)
      type is(fit_context_t)
         select case(ctx%family)
         case(FIT_HYPERPO);value=-loglik_hyperpo(theta(1:2),ctx%y)
         case(FIT_HYPERPO2);value=-loglik_hyperpo2(theta(1:2),ctx%y)
         case(FIT_DBH);value=loglik_dbh(logistic(theta(1)),ctx%y)
         case(FIT_DLD);value=loglik_dld(exp(theta(1)),ctx%y)
         case(FIT_DIKUM);value=-loglik_dikum(theta(1:2),ctx%y)
         case(FIT_POISXL);value=-loglik_poisxl(theta(1),ctx%y)
         case(FIT_GGEO);value=-loglik_ggeo(theta(1:2),ctx%y)
         case(FIT_DGEII);value=-loglik_dgeii(theta(1:2),ctx%y)
         case(FIT_DPERKS);value=-loglik_dperks(theta(1:2),ctx%y)
         case(FIT_DSPA);value=-loglik_dspa(theta(1:2),ctx%y)
         case(FIT_DMOLBE);value=-loglik_dmolbe(theta(1:2),ctx%y)
         case default;value=huge(1.0_dp)
         end select
      class default
         value=huge(1.0_dp)
      end select
      if(.not.(value<huge(1.0_dp)))value=huge(1.0_dp)/10.0_dp
   end function fit_objective

   subroutine run_fit(y,family,start,best,status)
      real(dp),intent(in)::y(:),start(:)
      integer,intent(in)::family
      real(dp),intent(out)::best(size(start))
      integer,intent(out),optional::status
      type(fit_context_t)::ctx
      real(dp)::fbest
      integer::st,it
      ctx%family=family;ctx%y=y
      call nelder_mead(fit_objective,ctx,start,best,fbest,st,it,max_iter=10000,tol=1.0e-10_dp,step=0.25_dp)
      if(present(status))status=st
   end subroutine run_fit

   function estim_mu_sigma_hyperpo(y,status) result(par)
      real(dp),intent(in)::y(:);integer,intent(out),optional::status
      real(dp)::par(2),best(2);integer::st
      call run_fit(y,FIT_HYPERPO,[0.0_dp,0.0_dp],best,st);par=exp(best);if(present(status))status=st
   end function estim_mu_sigma_hyperpo

   function estim_mu_sigma_hyperpo2(y,status) result(par)
      real(dp),intent(in)::y(:);integer,intent(out),optional::status
      real(dp)::par(2),best(2);integer::st
      call run_fit(y,FIT_HYPERPO2,[0.0_dp,0.0_dp],best,st);par=exp(best);if(present(status))status=st
   end function estim_mu_sigma_hyperpo2

   function estim_mu_dbh(y,status) result(mu)
      real(dp),intent(in)::y(:);integer,intent(out),optional::status
      real(dp)::mu,best(1);integer::st
      call run_fit(y,FIT_DBH,[0.0_dp],best,st);mu=logistic(best(1));if(present(status))status=st
   end function estim_mu_dbh

   function estim_mu_dld(y,status) result(mu)
      real(dp),intent(in)::y(:);integer,intent(out),optional::status
      real(dp)::mu,best(1);integer::st
      call run_fit(y,FIT_DLD,[log(0.5_dp)],best,st);mu=exp(best(1));if(present(status))status=st
   end function estim_mu_dld

   function estim_mu_sigma_dikum(y,status) result(par)
      real(dp),intent(in)::y(:);integer,intent(out),optional::status
      real(dp)::par(2),best(2);integer::st
      call run_fit(y,FIT_DIKUM,[0.0_dp,0.0_dp],best,st);par=exp(best);if(present(status))status=st
   end function estim_mu_sigma_dikum

   function estim_mu_poisxl(y,status) result(mu)
      real(dp),intent(in)::y(:);integer,intent(out),optional::status
      real(dp)::mu,best(1);integer::st
      call run_fit(y,FIT_POISXL,[0.0_dp],best,st);mu=exp(best(1));if(present(status))status=st
   end function estim_mu_poisxl

   function estim_mu_sigma_ggeo(y,status) result(par)
      real(dp),intent(in)::y(:);integer,intent(out),optional::status
      real(dp)::par(2),best(2);integer::st
      call run_fit(y,FIT_GGEO,[0.0_dp,0.0_dp],best,st);par=[logistic(best(1)),exp(best(2))]
      if(present(status))status=st
   end function estim_mu_sigma_ggeo

   function estim_mu_sigma_dgeii(y,status) result(par)
      real(dp),intent(in)::y(:);integer,intent(out),optional::status
      real(dp)::par(2),best(2),m,start1;integer::st
      m=max(0.0_dp,mean_sample(y));start1=1.0_dp-1.0_dp/(1.0_dp+m)
      call run_fit(y,FIT_DGEII,[start1,0.0_dp],best,st);par=[logistic(best(1)),exp(best(2))]
      if(present(status))status=st
   end function estim_mu_sigma_dgeii

   function estim_mu_sigma_dperks(y,status) result(par)
      real(dp),intent(in)::y(:);integer,intent(out),optional::status
      real(dp)::par(2),best(2);integer::st
      call run_fit(y,FIT_DPERKS,[0.0_dp,0.0_dp],best,st);par=exp(best);if(present(status))status=st
   end function estim_mu_sigma_dperks

   function estim_mu_sigma_dspa(y,status) result(par)
      real(dp),intent(in)::y(:);integer,intent(out),optional::status
      real(dp)::par(2),best(2);integer::st
      call run_fit(y,FIT_DSPA,[0.0_dp,0.0_dp],best,st);par=[exp(best(1)),logistic(best(2))]
      if(present(status))status=st
   end function estim_mu_sigma_dspa

   function estim_mu_sigma_dmolbe(y,status) result(par)
      real(dp),intent(in)::y(:);integer,intent(out),optional::status
      real(dp)::par(2),best(2);integer::st
      call run_fit(y,FIT_DMOLBE,[0.0_dp,0.0_dp],best,st);par=exp(best);if(present(status))status=st
   end function estim_mu_sigma_dmolbe

   function estim_mu_sigma_compo(y) result(par)
      real(dp),intent(in)::y(:)
      real(dp)::par(2)
      integer::imin,imax,nrange,i,k,best_start,best_len,cur_start,cur_len,npts
      integer,allocatable::counts(:)
      real(dp),allocatable::xx(:),yy(:)
      real(dp)::sx,sy,sxx,sxy,den,slope,intercept
      if(size(y)<2)then;par=[1.0_dp,1.0_dp];return;end if
      imin=int(minval(y));imax=int(maxval(y));nrange=imax-imin+1
      if(nrange<=1.or.nrange>1000000)then;par=[max(1.0e-6_dp,mean_sample(y)),1.0_dp];return;end if
      allocate(counts(nrange));counts=0
      do i=1,size(y);k=int(nint(y(i)))-imin+1;if(k>=1.and.k<=nrange)counts(k)=counts(k)+1;end do
      best_start=1;best_len=0;cur_start=1;cur_len=0
      do i=1,nrange
         if(counts(i)>0)then
            if(cur_len==0)cur_start=i
            cur_len=cur_len+1
            if(cur_len>best_len)then;best_len=cur_len;best_start=cur_start;end if
         else;cur_len=0;end if
      end do
      if(best_len<2)then;par=[max(1.0e-6_dp,mean_sample(y)),1.0_dp];return;end if
      npts=best_len-1;allocate(xx(npts),yy(npts))
      do i=1,npts
         if(imin+best_start+i-1<=0)then
            par=[max(1.0e-6_dp,mean_sample(y)),1.0_dp];return
         end if
         xx(i)=log(real(imin+best_start+i-1,dp))
         yy(i)=log(real(counts(best_start+i-1),dp)/real(counts(best_start+i),dp))
      end do
      sx=sum(xx);sy=sum(yy);sxx=sum(xx*xx);sxy=sum(xx*yy);den=real(npts,dp)*sxx-sx*sx
      if(abs(den)<1.0e-14_dp)then;par=[max(1.0e-6_dp,mean_sample(y)),1.0_dp];return;end if
      slope=(real(npts,dp)*sxy-sx*sy)/den;intercept=(sy-slope*sx)/real(npts,dp)
      par=[exp(-intercept),max(1.0e-6_dp,slope)]
   end function estim_mu_sigma_compo

end module discretedists_estimators
