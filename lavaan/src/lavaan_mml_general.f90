module lavaan_mml_general
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use lavaan_kinds, only : dp
   use lavaan_mml, only : gauss_hermite_normal
   use lavaan_ordinal, only : ordinal_thresholds
   use lavaan_linalg, only : chol_lower, inverse_general
   use lavaan_optimizer, only : bfgs_minimize
   use numderiv, only : hessian, nd_success
   implicit none
   private
   type, public :: mml_mixed_result
      real(dp),allocatable :: loadings(:,:), intercept(:), residual_sd(:), thresholds(:,:)
      real(dp),allocatable :: latent_mean(:), latent_cov(:,:), par(:), se(:), vcov(:,:)
      integer,allocatable :: ncat(:)
      logical,allocatable :: ordinal(:)
      real(dp) :: loglik=-huge(1.0_dp), aic=huge(1.0_dp), bic=huge(1.0_dp)
      integer :: status=0, iterations=0, nquad=0, nfactor=0
      logical :: converged=.false.
   end type mml_mixed_result
   public :: fit_mml_mixed_factor, mml_mixed_loglik
contains
   subroutine fit_mml_mixed_factor(data,ordinal,loading_start,free_mask,latent_mean,latent_cov,result,nquad)
      real(dp),intent(in)::data(:,:),loading_start(:,:),latent_mean(:),latent_cov(:,:)
      logical,intent(in)::ordinal(:),free_mask(:,:)
      type(mml_mixed_result),intent(out)::result
      integer,intent(in),optional::nquad
      real(dp),allocatable::x(:),hess(:,:),hi(:,:)
      integer::p,q,nq,kload,ncont,k,idx,i,j,info,status
      real(dp)::fval
      p=size(data,2)
      q=size(loading_start,2)
      nq=7
      if(present(nquad)) nq=max(3,nquad)
      if(size(ordinal)/=p .or. size(loading_start,1)/=p .or. any(shape(free_mask)/=shape(loading_start)) .or. &
         size(latent_mean)/=q .or. any(shape(latent_cov)/=[q,q]) .or. q<1 .or. q>3) then
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
      ! The likelihood combines quadrature and finite-difference gradients, so
      ! demanding substantially less than 1e-5 is below its practical noise floor.
      call bfgs_minimize(nll,x,fval,result%converged,result%iterations,maxiter=1600,tol=1.0e-5_dp)
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
      call hessian(nll,x,hess,status=status)
      allocate(result%vcov(k,k),result%se(k))
      result%vcov=0.0_dp
      result%se=huge(1.0_dp)
      if(status==nd_success) then
         call inverse_general(hess,hi,info)
         if(info==0) then
         result%vcov=hi
            do i=1,k
            if(hi(i,i)>=0.0_dp) result%se(i)=sqrt(hi(i,i))
            end do
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
         v=-mml_mixed_loglik(data,ordinal,ll,ii,rr,result%thresholds,result%ncat,latent_mean,latent_cov,nq)
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
   end subroutine fit_mml_mixed_factor

   function mml_mixed_loglik(data,ordinal,loadings,intercept,residual_sd,thresholds,ncat, &
                             latent_mean,latent_cov,nquad) result(ll)
      real(dp),intent(in)::data(:,:),loadings(:,:),intercept(:),residual_sd(:),thresholds(:,:)
      real(dp),intent(in)::latent_mean(:),latent_cov(:,:)
      logical,intent(in)::ordinal(:)
      integer,intent(in)::ncat(:),nquad
      real(dp)::ll
      real(dp),allocatable::nodes(:),weights(:),eta(:),z(:),l(:,:)
      integer::n,p,q,r,total,code,j,c,iq,tmp,info
      real(dp)::pr,cond,wt,lo,hi,mu,sdv,y,pi2
      n=size(data,1)
      p=size(data,2)
      q=size(loadings,2)
      ll=0.0_dp
      pi2=2.0_dp*acos(-1.0_dp)
      if(size(loadings,1)/=p .or. size(ordinal)/=p .or. size(intercept)/=p .or. size(residual_sd)/=p .or. &
         size(ncat)/=p .or. size(latent_mean)/=q .or. any(shape(latent_cov)/=[q,q]) .or. q<1 .or. q>3) then
         ll=-huge(1.0_dp)
         return
      end if
      call gauss_hermite_normal(nquad,nodes,weights,info)
      if(info/=0) then
      ll=-huge(1.0_dp)
      return
      end if
      call chol_lower(latent_cov,l,info)
      if(info/=0) then
      ll=-huge(1.0_dp)
      return
      end if
      allocate(eta(q),z(q))
      total=nquad**q
      do r=1,n
         pr=0.0_dp
         do code=0,total-1
            tmp=code
            wt=1.0_dp
            do iq=1,q
            c=mod(tmp,nquad)+1
            tmp=tmp/nquad
            z(iq)=nodes(c)
            wt=wt*weights(c)
            end do
            eta=latent_mean+matmul(l,z)
            cond=1.0_dp
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
                  cond=cond*max(normal_cdf(hi)-normal_cdf(lo),1.0e-300_dp)
               else
                  sdv=max(residual_sd(j),1.0e-8_dp)
                  y=data(r,j)
                  cond=cond*exp(-0.5_dp*((y-mu)/sdv)**2)/(sqrt(pi2)*sdv)
               end if
            end do
            pr=pr+wt*cond
         end do
         ll=ll+log(max(pr,1.0e-300_dp))
      end do
   end function mml_mixed_loglik

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
end module lavaan_mml_general
