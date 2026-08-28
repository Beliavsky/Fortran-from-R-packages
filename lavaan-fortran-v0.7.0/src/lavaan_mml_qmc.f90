module lavaan_mml_qmc
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use lavaan_kinds, only : dp
   use lavaan_linalg, only : chol_lower, inverse_general
   use lavaan_ordinal, only : ordinal_thresholds, normal_quantile
   use lavaan_optimizer, only : bfgs_minimize
   use lavaan_mml_general, only : mml_mixed_result
   use numderiv, only : hessian, nd_success
   implicit none
   private

   public :: mml_mixed_loglik_qmc, fit_mml_mixed_factor_qmc, halton_normal_nodes

contains

   function mml_mixed_loglik_qmc(data,ordinal,loadings,intercept,residual_sd,thresholds,ncat, &
                                 latent_mean,latent_cov,nsim) result(ll)
      real(dp),intent(in)::data(:,:),loadings(:,:),intercept(:),residual_sd(:),thresholds(:,:)
      real(dp),intent(in)::latent_mean(:),latent_cov(:,:)
      logical,intent(in)::ordinal(:)
      integer,intent(in)::ncat(:),nsim
      real(dp)::ll
      real(dp),allocatable::nodes(:,:),eta(:),l(:,:),logw(:)
      real(dp)::logcond,lo,hi,mu,sdv,y,pi2,mx,prob
      integer::n,p,q,r,k,j,c,info

      n=size(data,1)
      p=size(data,2)
      q=size(loadings,2)
      ll=0.0_dp
      pi2=2.0_dp*acos(-1.0_dp)
      if(size(loadings,1)/=p .or. size(ordinal)/=p .or. size(intercept)/=p .or. size(residual_sd)/=p .or. &
         size(ncat)/=p .or. size(latent_mean)/=q .or. any(shape(latent_cov)/=[q,q]) .or. &
         q<1 .or. q>8 .or. nsim<32) then
         ll=-huge(1.0_dp)
         return
      end if
      call chol_lower(latent_cov,l,info)
      if(info/=0) then
      ll=-huge(1.0_dp)
      return
      end if
      call halton_normal_nodes(q,nsim,nodes)
      allocate(eta(q),logw(nsim))
      do r=1,n
         do k=1,nsim
            eta=latent_mean+matmul(l,nodes(:,k))
            logcond=0.0_dp
            do j=1,p
               if(ieee_is_nan(data(r,j))) cycle
               mu=intercept(j)+dot_product(loadings(j,:),eta)
               if(ordinal(j)) then
                  c=nint(data(r,j))
                  if(c<1 .or. c>ncat(j)) then
                  ll=-huge(1.0_dp)
                  return
                  end if
                  if(c==1) then
                  lo=-huge(1.0_dp)
                  else
                  lo=thresholds(c-1,j)-mu
                  end if
                  if(c==ncat(j)) then
                  hi=huge(1.0_dp)
                  else
                  hi=thresholds(c,j)-mu
                  end if
                  prob=max(normal_cdf(hi)-normal_cdf(lo),1.0e-300_dp)
                  logcond=logcond+log(prob)
               else
                  sdv=max(residual_sd(j),1.0e-8_dp)
                  y=data(r,j)
                  logcond=logcond-0.5_dp*((y-mu)/sdv)**2-log(sqrt(pi2)*sdv)
               end if
            end do
            logw(k)=logcond
         end do
         mx=maxval(logw)
         ll=ll+mx+log(sum(exp(logw-mx))/real(nsim,dp))
      end do
   end function mml_mixed_loglik_qmc

   subroutine fit_mml_mixed_factor_qmc(data,ordinal,loading_start,free_mask,latent_mean,latent_cov,result,nsim,compute_se)
      real(dp),intent(in)::data(:,:),loading_start(:,:),latent_mean(:),latent_cov(:,:)
      logical,intent(in)::ordinal(:),free_mask(:,:)
      type(mml_mixed_result),intent(out)::result
      integer,intent(in),optional::nsim
      logical,intent(in),optional::compute_se
      real(dp),allocatable::x(:),hess(:,:),hi(:,:)
      integer::p,q,nq,kload,ncont,k,idx,i,j,info,status
      real(dp)::fval
      logical::dose

      p=size(data,2)
      q=size(loading_start,2)
      nq=1024
      if(present(nsim)) nq=max(64,nsim)
      dose=.true.
      if(present(compute_se)) dose=compute_se
      if(size(ordinal)/=p .or. size(loading_start,1)/=p .or. any(shape(free_mask)/=shape(loading_start)) .or. &
         size(latent_mean)/=q .or. any(shape(latent_cov)/=[q,q]) .or. q<1 .or. q>8) then
         result%status=-1
         return
      end if
      kload=count(free_mask)
      ncont=count(.not.ordinal)
      k=kload+2*ncont
      allocate(x(k))
      idx=0
      do j=1,q
         do i=1,p
            if(free_mask(i,j)) then
            idx=idx+1
            x(idx)=loading_start(i,j)
            end if
         end do
      end do
      allocate(result%intercept(p),result%residual_sd(p))
      result%intercept=0.0_dp
      result%residual_sd=1.0_dp
      do i=1,p
         if(.not.ordinal(i)) then
            idx=idx+1
            result%intercept(i)=observed_mean(data(:,i))
            x(idx)=result%intercept(i)
            idx=idx+1
            result%residual_sd(i)=max(0.1_dp,observed_sd(data(:,i),result%intercept(i)))
            x(idx)=log(result%residual_sd(i))
         end if
      end do
      call build_mixed_thresholds(data,ordinal,result%ncat,result%thresholds,info)
      if(info/=0) then
      result%status=info
      return
      end if
      call bfgs_minimize(nll,x,fval,result%converged,result%iterations,maxiter=1400,tol=3.0e-7_dp)
      result%loadings=loading_start
      call unpack(x,result%loadings,result%intercept,result%residual_sd)
      result%latent_mean=latent_mean
      result%latent_cov=latent_cov
      result%ordinal=ordinal
      result%par=x
      result%loglik=-fval
      result%aic=2.0_dp*fval+2.0_dp*real(k,dp)
      result%bic=2.0_dp*fval+log(real(size(data,1),dp))*real(k,dp)
      result%nquad=nq
      result%nfactor=q
      allocate(result%vcov(k,k),result%se(k))
      result%vcov=0.0_dp
      result%se=huge(1.0_dp)
      if(dose) then
         call hessian(nll,x,hess,status=status)
         if(status==nd_success) then
            call inverse_general(hess,hi,info)
            if(info==0) then
               result%vcov=hi
               do i=1,k
               if(hi(i,i)>=0.0_dp) result%se(i)=sqrt(hi(i,i))
               end do
            end if
         end if
      end if
      result%status=0
   contains
      function nll(z) result(v)
         real(dp),intent(in)::z(:)
         real(dp)::v
         real(dp),allocatable::ll(:,:),ii(:),rr(:)
         ll=loading_start
         ii=result%intercept
         rr=result%residual_sd
         call unpack(z,ll,ii,rr)
         v=-mml_mixed_loglik_qmc(data,ordinal,ll,ii,rr,result%thresholds,result%ncat,latent_mean,latent_cov,nq)
         if(.not.(v<huge(1.0_dp)/10.0_dp)) v=huge(1.0_dp)/100.0_dp
      end function nll
      subroutine unpack(z,ll,ii,rr)
         real(dp),intent(in)::z(:)
         real(dp),intent(inout)::ll(:,:),ii(:),rr(:)
         integer::a,b,pos
         pos=0
         do b=1,q
            do a=1,p
               if(free_mask(a,b)) then
               pos=pos+1
               ll(a,b)=z(pos)
               end if
            end do
         end do
         do a=1,p
            if(.not.ordinal(a)) then
               pos=pos+1
               ii(a)=z(pos)
               pos=pos+1
               rr(a)=exp(z(pos))
            end if
         end do
      end subroutine unpack
   end subroutine fit_mml_mixed_factor_qmc

   subroutine halton_normal_nodes(q,n,z)
      integer,intent(in)::q,n
      real(dp),allocatable,intent(out)::z(:,:)
      integer,parameter::primes(8)=[2,3,5,7,11,13,17,19]
      integer::i,j
      real(dp)::u
      allocate(z(q,n))
      do j=1,n
         do i=1,q
            u=radical_inverse(j+31,primes(i))
            u=min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,u))
            z(i,j)=normal_quantile(u)
         end do
      end do
   end subroutine halton_normal_nodes

   pure function radical_inverse(index,base) result(v)
      integer,intent(in)::index,base
      real(dp)::v,f
      integer::i
      i=index
      v=0.0_dp
      f=1.0_dp/real(base,dp)
      do while(i>0)
         v=v+f*real(mod(i,base),dp)
         i=i/base
         f=f/real(base,dp)
      end do
   end function radical_inverse

   subroutine build_mixed_thresholds(data,ordinal,ncat,thresholds,info)
      real(dp),intent(in)::data(:,:)
      logical,intent(in)::ordinal(:)
      integer,allocatable,intent(out)::ncat(:)
      real(dp),allocatable,intent(out)::thresholds(:,:)
      integer,intent(out)::info
      integer::p,j,r,c,maxc,nobs
      integer,allocatable::cnt(:)
      real(dp),allocatable::th(:)
      p=size(data,2)
      allocate(ncat(p))
      ncat=0
      maxc=1
      info=0
      do j=1,p
         if(ordinal(j)) then
            ncat(j)=0
            do r=1,size(data,1)
               if(.not.ieee_is_nan(data(r,j))) ncat(j)=max(ncat(j),nint(data(r,j)))
            end do
            if(ncat(j)<2) then
            info=j
            return
            end if
            maxc=max(maxc,ncat(j))
         end if
      end do
      allocate(thresholds(maxc-1,p))
      thresholds=0.0_dp
      do j=1,p
         if(.not.ordinal(j)) cycle
         allocate(cnt(ncat(j)))
         cnt=0
         nobs=0
         do r=1,size(data,1)
            if(ieee_is_nan(data(r,j))) cycle
            c=nint(data(r,j))
            if(c<1 .or. c>ncat(j)) then
            info=100+j
            return
            end if
            cnt(c)=cnt(c)+1
            nobs=nobs+1
         end do
         if(nobs==0 .or. any(cnt==0)) then
         info=200+j
         return
         end if
         th=ordinal_thresholds(cnt)
         thresholds(1:size(th),j)=th
         deallocate(cnt)
      end do
   end subroutine build_mixed_thresholds

   function observed_mean(x) result(m)
      real(dp),intent(in)::x(:)
      real(dp)::m
      integer::i,n
      m=0.0_dp
      n=0
      do i=1,size(x)
      if(.not.ieee_is_nan(x(i))) then
      m=m+x(i)
      n=n+1
      end if
      end do
      if(n>0) m=m/real(n,dp)
   end function observed_mean

   function observed_sd(x,m) result(s)
      real(dp),intent(in)::x(:),m
      real(dp)::s
      integer::i,n
      s=0.0_dp
      n=0
      do i=1,size(x)
      if(.not.ieee_is_nan(x(i))) then
      s=s+(x(i)-m)**2
      n=n+1
      end if
      end do
      if(n>1) then
      s=sqrt(s/real(n-1,dp))
      else
      s=1.0_dp
      end if
   end function observed_sd

   pure function normal_cdf(z) result(p)
      real(dp),intent(in)::z
      real(dp)::p
      if(z>8.0_dp) then
      p=1.0_dp
      else if(z< -8.0_dp) then
      p=0.0_dp
      else
      p=0.5_dp*erfc(-z/sqrt(2.0_dp))
      end if
   end function normal_cdf

end module lavaan_mml_qmc
