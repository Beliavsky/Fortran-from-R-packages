! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_glmrob
   use robustbase_kinds, only: dp
   use robustbase_linalg, only: solve_linear, invert_symmetric
   use robustbase_psi, only: huber_psi, huber_weight, tukey_psi, tukey_weight
   use robustbase_regression, only: robust_regression_result, robust_glm_fit
   implicit none
   private
   public :: glmrob_result, glmrob_mqle_fit, glmrob_mt_fit

   type :: glmrob_result
      real(dp),allocatable :: coefficients(:)
      real(dp),allocatable :: fitted(:)
      real(dp),allocatable :: residuals(:)
      real(dp),allocatable :: pearson_residuals(:)
      real(dp),allocatable :: weights(:)
      real(dp),allocatable :: covariance(:,:)
      real(dp),allocatable :: standard_errors(:)
      real(dp) :: objective=0.0_dp
      integer :: iterations=0
      character(len=8) :: method='Mqle'
      character(len=12) :: family='binomial'
      logical :: converged=.false.
   end type glmrob_result
contains
   subroutine glmrob_mqle_fit(x,y,family,result,tuning,start,max_iter,tol)
      real(dp),intent(in)::x(:,:),y(:)
      character(len=*),intent(in)::family
      type(glmrob_result),intent(out)::result
      real(dp),intent(in),optional::tuning,start(:),tol
      integer,intent(in),optional::max_iter
      type(robust_regression_result)::init
      integer::n,p,it,mi,info,i,half
      real(dp)::c,tt,delta,step,obj,newobj,q,epsi,edpsi,psi,score_i
      real(dp),allocatable::beta(:),candidate(:),eta(:),mu(:),var(:),dmu(:),r(:),weights(:),score(:),amat(:,:),bmat(:,:),ainv(:,:),xx(:,:),cov(:,:)
      logical::accepted
      n=size(x,1);p=size(x,2)
      if(size(y)/=n .or. p<1)error stop 'glmrob_mqle_fit: invalid dimensions'
      if(trim(family)=='binomial')then
         if(any(y<0.0_dp) .or. any(y>1.0_dp))error stop 'glmrob_mqle_fit: binomial y must be in [0,1]'
      else if(trim(family)=='poisson')then
         if(any(y<0.0_dp))error stop 'glmrob_mqle_fit: poisson y must be nonnegative'
      else
         error stop 'glmrob_mqle_fit: family must be binomial or poisson'
      end if
      c=1.345_dp
      if(present(tuning))c=tuning
      mi=100
      if(present(max_iter))mi=max(1,max_iter)
      tt=1.0e-7_dp
      if(present(tol))tt=max(tol,epsilon(1.0_dp))
      allocate(beta(p),candidate(p),eta(n),mu(n),var(n),dmu(n),r(n),weights(n),score(p),amat(p,p),bmat(p,p),ainv(p,p),xx(p,p),cov(p,p))
      if(present(start))then
         if(size(start)/=p)error stop 'glmrob_mqle_fit: start size mismatch'
         beta=start
      else
         call robust_glm_fit(x,y,family,init,max_iter=50)
         beta=init%coefficients
      end if
      call mean_variance(beta)
      obj=robust_objective(y,mu,var,c)
      result%converged=.false.
      do it=1,mi
         score=0.0_dp
         amat=0.0_dp
         do i=1,n
            call expected_huber(family,mu(i),var(i),c,epsi,edpsi)
            psi=huber_psi(r(i),c)
            q=dmu(i)/sqrt(var(i))
            score_i=q*(psi-epsi)
            score=score+x(i,:)*score_i
            xx=outer_product(x(i,:),x(i,:))
            amat=amat+q*q*max(edpsi,1.0e-8_dp)*xx
         end do
         call solve_linear(amat,score,candidate,info)
         if(info/=0)exit
         delta=maxval(abs(candidate))
         step=1.0_dp
         accepted=.false.
         do half=0,20
            candidate=beta+step*candidate
            call mean_variance_candidate(candidate,eta,mu,var,dmu,r)
            newobj=robust_objective(y,mu,var,c)
            if(newobj<=obj+1.0e-10_dp .or. step<=1.0e-5_dp)then
               accepted=.true.
               exit
            end if
            candidate=(candidate-beta)/step
            step=0.5_dp*step
         end do
         if(.not.accepted)exit
         beta=candidate
         obj=newobj
         if(delta*step<=tt*(1.0_dp+maxval(abs(beta))))then
            result%converged=.true.
            exit
         end if
      end do
      call mean_variance(beta)
      amat=0.0_dp;bmat=0.0_dp
      do i=1,n
         call expected_huber(family,mu(i),var(i),c,epsi,edpsi)
         psi=huber_psi(r(i),c)
         q=dmu(i)/sqrt(var(i))
         score_i=q*(psi-epsi)
         xx=outer_product(x(i,:),x(i,:))
         amat=amat+q*q*max(edpsi,1.0e-8_dp)*xx
         bmat=bmat+score_i*score_i*xx
         weights(i)=huber_weight(r(i),c)
      end do
      amat=amat/real(n,dp);bmat=bmat/real(n,dp)
      call invert_symmetric(amat,ainv,info,ridge=1.0e-10_dp)
      cov=matmul(matmul(ainv,bmat),ainv)/real(n,dp)
      cov=0.5_dp*(cov+transpose(cov))
      allocate(result%coefficients(p),result%fitted(n),result%residuals(n),result%pearson_residuals(n),result%weights(n),result%covariance(p,p),result%standard_errors(p))
      result%coefficients=beta;result%fitted=mu;result%residuals=y-mu;result%pearson_residuals=r;result%weights=weights;result%covariance=cov
      do i=1,p
         result%standard_errors(i)=sqrt(max(cov(i,i),0.0_dp))
      end do
      result%objective=obj;result%iterations=min(it,mi);result%method='Mqle';result%family=trim(family)
      if(info/=0)result%converged=.false.
   contains
      subroutine mean_variance(b)
         real(dp),intent(in)::b(:)
         call mean_variance_candidate(b,eta,mu,var,dmu,r)
      end subroutine mean_variance
      subroutine mean_variance_candidate(b,e,m,v,dm,rr)
         real(dp),intent(in)::b(:)
         real(dp),intent(out)::e(:),m(:),v(:),dm(:),rr(:)
         e=matmul(x,b)
         select case(trim(family))
         case('binomial')
            m=1.0_dp/(1.0_dp+exp(-max(-35.0_dp,min(35.0_dp,e))))
            v=max(m*(1.0_dp-m),1.0e-10_dp);dm=v
         case('poisson')
            m=exp(max(-25.0_dp,min(25.0_dp,e)))
            v=max(m,1.0e-10_dp);dm=m
         end select
         rr=(y-m)/sqrt(v)
      end subroutine mean_variance_candidate
   end subroutine glmrob_mqle_fit

   subroutine glmrob_mt_fit(x,y,result,trials,tuning,start,max_iter,tol,psi)
      real(dp),intent(in)::x(:,:),y(:)
      type(glmrob_result),intent(out)::result
      real(dp),intent(in),optional::trials(:),tuning,start(:),tol
      integer,intent(in),optional::max_iter
      character(len=*),intent(in),optional::psi
      type(robust_regression_result)::init
      integer::n,p,it,mi,info,i
      real(dp)::c,tt,delta
      character(len=8)::ps
      real(dp),allocatable::m(:),z(:),beta(:),newbeta(:),eta(:),mu(:),tmu(:),res(:),jac(:,:),w(:),a(:,:),bmat(:,:),ainv(:,:),xx(:,:),cov(:,:)
      n=size(x,1);p=size(x,2)
      if(size(y)/=n .or. p<1)error stop 'glmrob_mt_fit: invalid dimensions'
      allocate(m(n),z(n),beta(p),newbeta(p),eta(n),mu(n),tmu(n),res(n),jac(n,p),w(n),a(p,p),bmat(p,p),ainv(p,p),xx(p,p),cov(p,p))
      m=1.0_dp
      if(present(trials))then
         if(size(trials)/=n)error stop 'glmrob_mt_fit: trials size mismatch'
         m=trials
      end if
      if(any(m<=0.0_dp) .or. any(y<0.0_dp) .or. any(y>m))error stop 'glmrob_mt_fit: invalid binomial counts'
      c=4.685061_dp
      if(present(tuning))c=tuning
      ps='tukey'
      if(present(psi))ps=adjustl(psi)
      mi=100
      if(present(max_iter))mi=max(1,max_iter)
      tt=1.0e-7_dp
      if(present(tol))tt=max(tol,epsilon(1.0_dp))
      z=2.0_dp*asin(sqrt((y+0.375_dp)/(m+0.75_dp)))
      if(present(start))then
         if(size(start)/=p)error stop 'glmrob_mt_fit: start size mismatch'
         beta=start
      else
         call robust_glm_fit(x,y/m,'binomial',init,max_iter=50)
         beta=init%coefficients
      end if
      info=0
      do it=1,mi
         eta=matmul(x,beta)
         mu=1.0_dp/(1.0_dp+exp(-max(-35.0_dp,min(35.0_dp,eta))))
         tmu=2.0_dp*asin(sqrt(mu))
         res=(z-tmu)*sqrt(m+0.5_dp)
         do i=1,p
            jac(:,i)=x(:,i)*sqrt(max(mu*(1.0_dp-mu),1.0e-12_dp))*sqrt(m+0.5_dp)
         end do
         select case(trim(ps))
         case('huber')
            w=huber_weight(res,c)
         case default
            w=tukey_weight(res,c)
         end select
         call weighted_step(jac,res,w,newbeta,info)
         if(info/=0)exit
         delta=maxval(abs(newbeta))
         beta=beta+newbeta
         if(delta<=tt*(1.0_dp+maxval(abs(beta))))then
            result%converged=.true.
            exit
         end if
      end do
      eta=matmul(x,beta);mu=1.0_dp/(1.0_dp+exp(-max(-35.0_dp,min(35.0_dp,eta))));tmu=2.0_dp*asin(sqrt(mu));res=(z-tmu)*sqrt(m+0.5_dp)
      do i=1,p
         jac(:,i)=x(:,i)*sqrt(max(mu*(1.0_dp-mu),1.0e-12_dp))*sqrt(m+0.5_dp)
      end do
      if(trim(ps)=='huber')then
         w=huber_weight(res,c)
      else
         w=tukey_weight(res,c)
      end if
      a=0.0_dp;bmat=0.0_dp
      do i=1,n
         xx=outer_product(jac(i,:),jac(i,:))
         if(trim(ps)=='huber')then
            a=a+merge(1.0_dp,0.0_dp,abs(res(i))<c)*xx
            bmat=bmat+huber_psi(res(i),c)**2*xx
         else
            a=a+tukey_derivative(res(i),c)*xx
            bmat=bmat+tukey_psi(res(i),c)**2*xx
         end if
      end do
      a=a/real(n,dp);bmat=bmat/real(n,dp)
      call invert_symmetric(a,ainv,info,ridge=1.0e-10_dp)
      cov=matmul(matmul(ainv,bmat),ainv)/real(n,dp);cov=0.5_dp*(cov+transpose(cov))
      allocate(result%coefficients(p),result%fitted(n),result%residuals(n),result%pearson_residuals(n),result%weights(n),result%covariance(p,p),result%standard_errors(p))
      result%coefficients=beta;result%fitted=m*mu;result%residuals=y-result%fitted;result%pearson_residuals=res;result%weights=w;result%covariance=cov
      do i=1,p
         result%standard_errors(i)=sqrt(max(cov(i,i),0.0_dp))
      end do
      result%objective=sum(merge(0.5_dp*res*res,c*abs(res)-0.5_dp*c*c,abs(res)<=c))
      result%iterations=min(it,mi);result%method='MT';result%family='binomial'
      if(info/=0)result%converged=.false.
   end subroutine glmrob_mt_fit

   subroutine expected_huber(family,mu,var,c,epsi,edpsi)
      character(len=*),intent(in)::family
      real(dp),intent(in)::mu,var,c
      real(dp),intent(out)::epsi,edpsi
      real(dp)::p,r,pk,total
      integer::k,kmax
      if(trim(family)=='binomial')then
         r=(0.0_dp-mu)/sqrt(var)
         epsi=(1.0_dp-mu)*huber_psi(r,c)
         edpsi=(1.0_dp-mu)*merge(1.0_dp,0.0_dp,abs(r)<c)
         r=(1.0_dp-mu)/sqrt(var)
         epsi=epsi+mu*huber_psi(r,c)
         edpsi=edpsi+mu*merge(1.0_dp,0.0_dp,abs(r)<c)
      else
         if(mu>80.0_dp)then
            epsi=0.0_dp
            edpsi=erf(c/sqrt(2.0_dp))
            return
         end if
         pk=exp(-mu);total=pk;epsi=0.0_dp;edpsi=0.0_dp
         kmax=max(30,int(mu+12.0_dp*sqrt(max(mu,1.0_dp))+20.0_dp))
         do k=0,kmax
            if(k>0)then
               pk=pk*mu/real(k,dp)
               total=total+pk
            end if
            r=(real(k,dp)-mu)/sqrt(var)
            epsi=epsi+pk*huber_psi(r,c)
            edpsi=edpsi+pk*merge(1.0_dp,0.0_dp,abs(r)<c)
         end do
         p=max(0.0_dp,1.0_dp-total)
         if(p>1.0e-12_dp)then
            r=(real(kmax+1,dp)-mu)/sqrt(var)
            epsi=epsi+p*huber_psi(r,c)
            edpsi=edpsi+p*merge(1.0_dp,0.0_dp,abs(r)<c)
         end if
      end if
   end subroutine expected_huber

   function robust_objective(y,mu,var,c) result(value)
      real(dp),intent(in)::y(:),mu(:),var(:),c
      real(dp)::value
      real(dp),allocatable::r(:)
      r=(y-mu)/sqrt(var)
      value=sum(merge(0.5_dp*r*r,c*abs(r)-0.5_dp*c*c,abs(r)<=c))
   end function robust_objective

   subroutine weighted_step(jac,res,w,step,info)
      real(dp),intent(in)::jac(:,:),res(:),w(:)
      real(dp),intent(out)::step(:)
      integer,intent(out)::info
      real(dp),allocatable::a(:,:),b(:)
      integer::j
      allocate(a(size(jac,1),size(jac,2)),b(size(res)))
      do j=1,size(jac,2)
         a(:,j)=jac(:,j)*sqrt(max(w,0.0_dp))
      end do
      b=res*sqrt(max(w,0.0_dp))
      call least_squares_local(a,b,step,info)
   end subroutine weighted_step

   subroutine least_squares_local(a,b,x,info)
      use robustbase_linalg, only: least_squares
      real(dp),intent(in)::a(:,:),b(:)
      real(dp),intent(out)::x(:)
      integer,intent(out)::info
      call least_squares(a,b,x,info)
   end subroutine least_squares_local

   elemental function tukey_derivative(x,c) result(value)
      real(dp),intent(in)::x,c
      real(dp)::value,u
      u=x/c
      if(abs(u)<1.0_dp)then
         value=(1.0_dp-u*u)*(1.0_dp-5.0_dp*u*u)
      else
         value=0.0_dp
      end if
   end function tukey_derivative

   pure function outer_product(a,b) result(c)
      real(dp),intent(in)::a(:),b(:)
      real(dp)::c(size(a),size(b))
      integer::j
      do j=1,size(b)
         c(:,j)=a*b(j)
      end do
   end function outer_product
end module robustbase_glmrob
