! Generic censored GAMLSS likelihoods using translated gamlss.dist CDFs.
! Matrix-first counterpart of cens()/Surv handling in upstream gamlss.
! v0.4 uses a module-level optimizer context so no nested-procedure
! trampoline or executable stack is required by GNU Fortran.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_censoring
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use gamlss_kinds, only : dp
   use gamlss_fit, only : gamlss_fit_result_t, family_npar, map_parameters, &
      family_logpdf, default_parameters, inverse_link
   use gamlss_optim, only : bfgs_minimize, numerical_hessian
   use gamlss_linalg, only : invert_matrix
   use gamlss_family_support, only : family_cdf, family_is_discrete
   implicit none
   private
   integer,parameter,public :: CENS_EXACT=0,CENS_LEFT=1,CENS_RIGHT=2,CENS_INTERVAL=3
   public :: fit_gamlss_censored, censored_case_loglik, truncated_censored_case_loglik
   public :: surv_right_censoring, surv_interval2, surv_counting_process

   type :: censor_context_t
      integer :: n=0,family=0,p1=0,p2=0,p3=0,p4=0,o2=0,o3=0,o4=0,ntheta=0
      real(dp),allocatable :: lower(:),upper(:),x_mu(:,:),x_sigma(:,:),x_nu(:,:),x_tau(:,:)
      real(dp),allocatable :: weights(:),off_mu(:),off_sigma(:),off_nu(:),off_tau(:),entry(:)
      integer,allocatable :: censor(:)
      logical :: has_entry=.false.
   end type censor_context_t
   type(censor_context_t),save :: censor_ctx
contains

   subroutine surv_right_censoring(time,event,lower,upper,censor,status)
      real(dp),intent(in)::time(:)
      integer,intent(in)::event(:)
      real(dp),allocatable,intent(out)::lower(:),upper(:)
      integer,allocatable,intent(out)::censor(:)
      integer,intent(out),optional::status
      integer::i,istat
      istat=0
      if(size(event)/=size(time).or.any(event<0).or.any(event>1))then
         allocate(lower(0),upper(0),censor(0));istat=1
         if(present(status))status=istat
         return
      end if
      allocate(lower(size(time)),upper(size(time)),censor(size(time)))
      do i=1,size(time)
         lower(i)=time(i);upper(i)=time(i)
         if(event(i)==1)then
            censor(i)=CENS_EXACT
         else
            censor(i)=CENS_RIGHT
         end if
      end do
      if(present(status))status=istat
   end subroutine surv_right_censoring

   subroutine surv_interval2(lower_in,upper_in,lower,upper,censor,status)
      real(dp),intent(in)::lower_in(:),upper_in(:)
      real(dp),allocatable,intent(out)::lower(:),upper(:)
      integer,allocatable,intent(out)::censor(:)
      integer,intent(out),optional::status
      integer::i,istat
      logical::lf,uf
      istat=0
      if(size(upper_in)/=size(lower_in))then
         allocate(lower(0),upper(0),censor(0));istat=1
         if(present(status))status=istat
         return
      end if
      allocate(lower(size(lower_in)),upper(size(lower_in)),censor(size(lower_in)))
      lower=lower_in;upper=upper_in
      do i=1,size(lower_in)
         lf=ieee_is_finite(lower_in(i));uf=ieee_is_finite(upper_in(i))
         if(lf.and.uf)then
            if(upper_in(i)<lower_in(i))then;istat=2;exit;end if
            if(abs(upper_in(i)-lower_in(i))<=1.0e-14_dp)then
               censor(i)=CENS_EXACT
            else
               censor(i)=CENS_INTERVAL
            end if
         else if(lf.and..not.uf)then
            censor(i)=CENS_RIGHT;upper(i)=lower(i)
         else if(.not.lf.and.uf)then
            censor(i)=CENS_LEFT;lower(i)=upper(i)
         else
            istat=3;exit
         end if
      end do
      if(present(status))status=istat
   end subroutine surv_interval2

   subroutine surv_counting_process(start,stop,event,entry,lower,upper,censor,status)
      real(dp),intent(in)::start(:),stop(:)
      integer,intent(in)::event(:)
      real(dp),allocatable,intent(out)::entry(:),lower(:),upper(:)
      integer,allocatable,intent(out)::censor(:)
      integer,intent(out),optional::status
      integer::i,istat
      istat=0
      if(size(stop)/=size(start).or.size(event)/=size(start).or.any(stop<=start) .or. &
         any(event<0).or.any(event>1))then
         allocate(entry(0),lower(0),upper(0),censor(0));istat=1
         if(present(status))status=istat
         return
      end if
      allocate(entry(size(start)),lower(size(start)),upper(size(start)),censor(size(start)))
      entry=start;lower=stop;upper=stop
      do i=1,size(start)
         censor(i)=merge(CENS_EXACT,CENS_RIGHT,event(i)==1)
      end do
      if(present(status))status=istat
   end subroutine surv_counting_process

   subroutine fit_gamlss_censored(lower,upper,censor,x_mu,family,result,x_sigma,x_nu,x_tau, &
      weights,offset_mu,offset_sigma,offset_nu,offset_tau,start,max_iter,tol,entry)
      real(dp),intent(in)::lower(:),upper(:),x_mu(:,:)
      integer,intent(in)::censor(:),family
      type(gamlss_fit_result_t),intent(out)::result
      real(dp),intent(in),optional::x_sigma(:,:),x_nu(:,:),x_tau(:,:),weights(:)
      real(dp),intent(in),optional::offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:),start(:),entry(:)
      integer,intent(in),optional::max_iter
      real(dp),intent(in),optional::tol
      integer::n,np,p1,p2,p3,p4,ntheta,o2,o3,o4,istat,i
      real(dp),allocatable::theta(:),hess(:,:),cov(:,:),w(:),pseudo(:)
      real(dp)::fval,tolerance

      n=size(lower); np=family_npar(family); p1=size(x_mu,2)
      if(n<=0 .or. size(upper)/=n .or. size(censor)/=n .or. size(x_mu,1)/=n .or. np<1)then
         result%status=10; return
      end if
      if(any(censor<CENS_EXACT).or.any(censor>CENS_INTERVAL))then;result%status=11;return;end if
      if(any(censor==CENS_INTERVAL .and. upper<lower))then;result%status=12;return;end if
      p2=0;p3=0;p4=0
      if(np>=2)then
         if(.not.present(x_sigma))then;result%status=13;return;end if
         if(size(x_sigma,1)/=n)then;result%status=14;return;end if
         p2=size(x_sigma,2)
      end if
      if(np>=3)then
         if(.not.present(x_nu))then;result%status=15;return;end if
         if(size(x_nu,1)/=n)then;result%status=16;return;end if
         p3=size(x_nu,2)
      end if
      if(np>=4)then
         if(.not.present(x_tau))then;result%status=17;return;end if
         if(size(x_tau,1)/=n)then;result%status=18;return;end if
         p4=size(x_tau,2)
      end if
      if(present(offset_mu))then;if(size(offset_mu)/=n)then;result%status=19;return;end if;end if
      if(present(offset_sigma))then;if(size(offset_sigma)/=n)then;result%status=20;return;end if;end if
      if(present(offset_nu))then;if(size(offset_nu)/=n)then;result%status=21;return;end if;end if
      if(present(offset_tau))then;if(size(offset_tau)/=n)then;result%status=22;return;end if;end if
      if(present(entry))then;if(size(entry)/=n)then;result%status=25;return;end if;end if
      allocate(w(n));w=1.0_dp
      if(present(weights))then
         if(size(weights)/=n.or.any(weights<0.0_dp))then;result%status=23;return;end if
         w=weights
      end if
      ntheta=p1+p2+p3+p4;o2=p1;o3=p1+p2;o4=p1+p2+p3
      allocate(theta(ntheta));theta=0.0_dp
      call set_censor_context(lower,upper,censor,x_mu,family,w,x_sigma,x_nu,x_tau, &
         offset_mu,offset_sigma,offset_nu,offset_tau,entry)
      if(present(start))then
         if(size(start)/=ntheta)then;result%status=24;call clear_censor_context();return;end if
         theta=start
      else
         allocate(pseudo(n))
         do i=1,n
            select case(censor(i))
            case(CENS_EXACT,CENS_RIGHT); pseudo(i)=lower(i)
            case(CENS_LEFT); pseudo(i)=upper(i)
            case(CENS_INTERVAL); pseudo(i)=0.5_dp*(lower(i)+upper(i))
            end select
         end do
         call initialize_censor_parameters(theta,pseudo)
      end if
      tolerance=1.0e-7_dp;if(present(tol))tolerance=tol
      call bfgs_minimize(censor_objective,theta,fval,result%iterations_status,max_iter=max_iter,tol=tolerance)
      result%converged=(result%iterations_status==0)
      result%status=result%iterations_status;result%loglik=-fval
      result%aic=2.0_dp*real(ntheta,dp)-2.0_dp*result%loglik
      result%beta_mu=theta(1:p1)
      if(p2>0)result%beta_sigma=theta(o2+1:o3)
      if(p3>0)result%beta_nu=theta(o3+1:o4)
      if(p4>0)result%beta_tau=theta(o4+1:ntheta)
      call fill_censor_fitted(theta,result)
      allocate(hess(ntheta,ntheta));call numerical_hessian(censor_objective,theta,hess)
      call invert_matrix(hess,cov,istat)
      if(istat==0)then;result%covariance=cov;else;allocate(result%covariance(0,0));end if
      call clear_censor_context()
   end subroutine fit_gamlss_censored

   subroutine set_censor_context(lower,upper,censor,x_mu,family,w,x_sigma,x_nu,x_tau, &
      offset_mu,offset_sigma,offset_nu,offset_tau,entry)
      real(dp),intent(in)::lower(:),upper(:),x_mu(:,:),w(:)
      integer,intent(in)::censor(:),family
      real(dp),intent(in),optional::x_sigma(:,:),x_nu(:,:),x_tau(:,:)
      real(dp),intent(in),optional::offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:),entry(:)
      integer::n
      call clear_censor_context()
      n=size(lower);censor_ctx%n=n;censor_ctx%family=family;censor_ctx%p1=size(x_mu,2)
      censor_ctx%p2=0;censor_ctx%p3=0;censor_ctx%p4=0
      if(present(x_sigma))censor_ctx%p2=size(x_sigma,2)
      if(present(x_nu))censor_ctx%p3=size(x_nu,2)
      if(present(x_tau))censor_ctx%p4=size(x_tau,2)
      censor_ctx%o2=censor_ctx%p1;censor_ctx%o3=censor_ctx%p1+censor_ctx%p2
      censor_ctx%o4=censor_ctx%o3+censor_ctx%p3
      censor_ctx%ntheta=censor_ctx%o4+censor_ctx%p4
      censor_ctx%lower=lower;censor_ctx%upper=upper;censor_ctx%censor=censor;censor_ctx%x_mu=x_mu
      censor_ctx%weights=w
      if(present(x_sigma))then;censor_ctx%x_sigma=x_sigma;else;allocate(censor_ctx%x_sigma(n,0));end if
      if(present(x_nu))then;censor_ctx%x_nu=x_nu;else;allocate(censor_ctx%x_nu(n,0));end if
      if(present(x_tau))then;censor_ctx%x_tau=x_tau;else;allocate(censor_ctx%x_tau(n,0));end if
      allocate(censor_ctx%off_mu(n),censor_ctx%off_sigma(n),censor_ctx%off_nu(n),censor_ctx%off_tau(n))
      censor_ctx%off_mu=0.0_dp;censor_ctx%off_sigma=0.0_dp;censor_ctx%off_nu=0.0_dp;censor_ctx%off_tau=0.0_dp
      if(present(offset_mu))censor_ctx%off_mu=offset_mu
      if(present(offset_sigma))censor_ctx%off_sigma=offset_sigma
      if(present(offset_nu))censor_ctx%off_nu=offset_nu
      if(present(offset_tau))censor_ctx%off_tau=offset_tau
      censor_ctx%has_entry=present(entry)
      if(present(entry))then;censor_ctx%entry=entry;else;allocate(censor_ctx%entry(0));end if
   end subroutine set_censor_context

   subroutine clear_censor_context()
      if(allocated(censor_ctx%lower))deallocate(censor_ctx%lower)
      if(allocated(censor_ctx%upper))deallocate(censor_ctx%upper)
      if(allocated(censor_ctx%censor))deallocate(censor_ctx%censor)
      if(allocated(censor_ctx%x_mu))deallocate(censor_ctx%x_mu)
      if(allocated(censor_ctx%x_sigma))deallocate(censor_ctx%x_sigma)
      if(allocated(censor_ctx%x_nu))deallocate(censor_ctx%x_nu)
      if(allocated(censor_ctx%x_tau))deallocate(censor_ctx%x_tau)
      if(allocated(censor_ctx%weights))deallocate(censor_ctx%weights)
      if(allocated(censor_ctx%off_mu))deallocate(censor_ctx%off_mu)
      if(allocated(censor_ctx%off_sigma))deallocate(censor_ctx%off_sigma)
      if(allocated(censor_ctx%off_nu))deallocate(censor_ctx%off_nu)
      if(allocated(censor_ctx%off_tau))deallocate(censor_ctx%off_tau)
      if(allocated(censor_ctx%entry))deallocate(censor_ctx%entry)
      censor_ctx%n=0;censor_ctx%family=0;censor_ctx%has_entry=.false.
   end subroutine clear_censor_context

   subroutine initialize_censor_parameters(par,yy)
      real(dp),intent(out)::par(:)
      real(dp),intent(in)::yy(:)
      real(dp)::m,s,a,b,c,d
      par=0.0_dp;m=sum(yy)/real(censor_ctx%n,dp)
      s=sqrt(max(sum((yy-m)**2)/real(max(1,censor_ctx%n-1),dp),1.0e-8_dp))
      call default_parameters(censor_ctx%family,m,s,a,b,c,d)
      if(is_intercept(censor_ctx%x_mu(:,1)))par(1)=inverse_link(censor_ctx%family,1,a)-mean_vec(censor_ctx%off_mu)
      if(censor_ctx%p2>0)then
         if(is_intercept(censor_ctx%x_sigma(:,1)))then
            par(censor_ctx%o2+1)=inverse_link(censor_ctx%family,2,b)-mean_vec(censor_ctx%off_sigma)
         end if
      end if
      if(censor_ctx%p3>0)then
         if(is_intercept(censor_ctx%x_nu(:,1)))then
            par(censor_ctx%o3+1)=inverse_link(censor_ctx%family,3,c)-mean_vec(censor_ctx%off_nu)
         end if
      end if
      if(censor_ctx%p4>0)then
         if(is_intercept(censor_ctx%x_tau(:,1)))then
            par(censor_ctx%o4+1)=inverse_link(censor_ctx%family,4,d)-mean_vec(censor_ctx%off_tau)
         end if
      end if
   end subroutine initialize_censor_parameters

   real(dp) function mean_vec(x) result(v)
      real(dp),intent(in)::x(:)
      v=sum(x)/real(max(1,size(x)),dp)
   end function mean_vec

   logical function is_intercept(col) result(ok)
      real(dp),intent(in)::col(:)
      ok=maxval(abs(col-1.0_dp))<1.0e-10_dp
   end function is_intercept

   real(dp) function censor_objective(par) result(nll)
      real(dp),intent(in)::par(:)
      real(dp)::e1,e2,e3,e4,a,b,c,d,ll
      integer::i
      nll=0.0_dp
      do i=1,censor_ctx%n
         e1=dot_product(censor_ctx%x_mu(i,:),par(1:censor_ctx%p1))+censor_ctx%off_mu(i)
         e2=0.0_dp;e3=0.0_dp;e4=0.0_dp
         if(censor_ctx%p2>0)then
            e2=dot_product(censor_ctx%x_sigma(i,:),par(censor_ctx%o2+1:censor_ctx%o3))+censor_ctx%off_sigma(i)
         end if
         if(censor_ctx%p3>0)then
            e3=dot_product(censor_ctx%x_nu(i,:),par(censor_ctx%o3+1:censor_ctx%o4))+censor_ctx%off_nu(i)
         end if
         if(censor_ctx%p4>0)then
            e4=dot_product(censor_ctx%x_tau(i,:),par(censor_ctx%o4+1:censor_ctx%ntheta))+censor_ctx%off_tau(i)
         end if
         call map_parameters(censor_ctx%family,e1,e2,e3,e4,a,b,c,d)
         if(censor_ctx%has_entry)then
            ll=truncated_censored_case_loglik(censor_ctx%family,censor_ctx%lower(i),censor_ctx%upper(i), &
               censor_ctx%censor(i),censor_ctx%entry(i),a,b,c,d)
         else
            ll=censored_case_loglik(censor_ctx%family,censor_ctx%lower(i),censor_ctx%upper(i), &
               censor_ctx%censor(i),a,b,c,d)
         end if
         if(.not.ieee_is_finite(ll))then
            nll=nll+1.0e12_dp+1.0e4_dp*sum(par*par)
         else
            nll=nll-censor_ctx%weights(i)*ll
         end if
      end do
   end function censor_objective

   subroutine fill_censor_fitted(par,result)
      real(dp),intent(in)::par(:)
      type(gamlss_fit_result_t),intent(inout)::result
      real(dp)::e1,e2,e3,e4,a,b,c,d
      integer::i
      allocate(result%fitted_mu(censor_ctx%n))
      if(censor_ctx%p2>0)allocate(result%fitted_sigma(censor_ctx%n))
      if(censor_ctx%p3>0)allocate(result%fitted_nu(censor_ctx%n))
      if(censor_ctx%p4>0)allocate(result%fitted_tau(censor_ctx%n))
      do i=1,censor_ctx%n
         e1=dot_product(censor_ctx%x_mu(i,:),par(1:censor_ctx%p1))+censor_ctx%off_mu(i)
         e2=0.0_dp;e3=0.0_dp;e4=0.0_dp
         if(censor_ctx%p2>0)then
            e2=dot_product(censor_ctx%x_sigma(i,:),par(censor_ctx%o2+1:censor_ctx%o3))+censor_ctx%off_sigma(i)
         end if
         if(censor_ctx%p3>0)then
            e3=dot_product(censor_ctx%x_nu(i,:),par(censor_ctx%o3+1:censor_ctx%o4))+censor_ctx%off_nu(i)
         end if
         if(censor_ctx%p4>0)then
            e4=dot_product(censor_ctx%x_tau(i,:),par(censor_ctx%o4+1:censor_ctx%ntheta))+censor_ctx%off_tau(i)
         end if
         call map_parameters(censor_ctx%family,e1,e2,e3,e4,a,b,c,d)
         result%fitted_mu(i)=a
         if(censor_ctx%p2>0)result%fitted_sigma(i)=b
         if(censor_ctx%p3>0)result%fitted_nu(i)=c
         if(censor_ctx%p4>0)result%fitted_tau(i)=d
      end do
   end subroutine fill_censor_fitted

   real(dp) function censored_case_loglik(family,lower,upper,code,a,b,c,d) result(ll)
      integer,intent(in)::family,code
      real(dp),intent(in)::lower,upper,a,b,c,d
      real(dp)::fl,fu,prob,ql,qu
      logical::disc
      disc=family_is_discrete(family)
      select case(code)
      case(CENS_EXACT)
         ll=family_logpdf(family,lower,a,b,c,d)
      case(CENS_LEFT)
         qu=upper;if(disc)qu=floor(upper)
         fu=family_cdf(family,qu,a,b,c,d)
         if(fu<0.0_dp)then;ll=-huge(1.0_dp);else;ll=log(max(tiny(1.0_dp),fu));end if
      case(CENS_RIGHT)
         ql=lower;if(disc)ql=floor(lower)
         fl=family_cdf(family,ql,a,b,c,d)
         if(fl<0.0_dp)then;ll=-huge(1.0_dp);else;ll=log(max(tiny(1.0_dp),1.0_dp-fl));end if
      case(CENS_INTERVAL)
         ql=lower;qu=upper
         if(disc)then;ql=floor(lower);qu=floor(upper);end if
         fl=family_cdf(family,ql,a,b,c,d);fu=family_cdf(family,qu,a,b,c,d)
         if(fl<0.0_dp.or.fu<0.0_dp)then
            ll=-huge(1.0_dp)
         else
            prob=max(tiny(1.0_dp),fu-fl);ll=log(prob)
         end if
      case default
         ll=-huge(1.0_dp)
      end select
   end function censored_case_loglik

   real(dp) function truncated_censored_case_loglik(family,lower,upper,code,entry,a,b,c,d) result(ll)
      integer,intent(in)::family,code
      real(dp),intent(in)::lower,upper,entry,a,b,c,d
      real(dp)::fe,fl,fu,den,prob,qe,ql,qu
      logical::disc
      disc=family_is_discrete(family)
      qe=entry;if(disc)qe=floor(entry)
      fe=family_cdf(family,qe,a,b,c,d)
      if(fe<0.0_dp)then;ll=-huge(1.0_dp);return;end if
      den=max(tiny(1.0_dp),1.0_dp-fe)
      select case(code)
      case(CENS_EXACT)
         if(lower<=entry)then;ll=-huge(1.0_dp);return;end if
         ll=family_logpdf(family,lower,a,b,c,d)-log(den)
      case(CENS_LEFT)
         qu=upper;if(disc)qu=floor(upper)
         if(upper<=entry)then;ll=-huge(1.0_dp);return;end if
         fu=family_cdf(family,qu,a,b,c,d)
         prob=max(tiny(1.0_dp),fu-fe);ll=log(prob)-log(den)
      case(CENS_RIGHT)
         ql=max(lower,entry);if(disc)ql=floor(ql)
         fl=family_cdf(family,ql,a,b,c,d)
         prob=max(tiny(1.0_dp),1.0_dp-fl);ll=log(prob)-log(den)
      case(CENS_INTERVAL)
         if(upper<=entry)then;ll=-huge(1.0_dp);return;end if
         ql=max(lower,entry);qu=upper
         if(disc)then;ql=floor(ql);qu=floor(qu);end if
         fl=family_cdf(family,ql,a,b,c,d);fu=family_cdf(family,qu,a,b,c,d)
         prob=max(tiny(1.0_dp),fu-fl);ll=log(prob)-log(den)
      case default
         ll=-huge(1.0_dp)
      end select
   end function truncated_censored_case_loglik

end module gamlss_censoring
