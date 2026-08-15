! Residual and model-diagnostic computations from gamlss.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_diagnostics
   use gamlss_kinds, only : dp
   use gamlss_fit
   use gamlss_continuous, only : pNO,pGA,pBE,pWEI,pLNO,pTF,pEGB2,pGB2,pGG,pJSUo
   use gamlss_discrete, only : pNBI,pNBII,pZIP,pPIG
   use gamlss_boxcox, only : pBCCG,pBCT,pBCPE
   use gamlss_special, only : normal_quantile
   use gamlss_types, only : gamlss_result_t
   use gamlss_linalg, only : invert_matrix
   implicit none
   private
   public :: randomized_quantile_residuals, residual_moments, residual_acf
   public :: gaic, likelihood_ratio_stat, deviance_increment
   public :: hat_values_penalized, effective_df_penalized

contains

   subroutine randomized_quantile_residuals(y,result,residuals,status,randomize)
      real(dp),intent(in)::y(:)
      type(gamlss_result_t),intent(in)::result
      real(dp),allocatable,intent(out)::residuals(:)
      integer,intent(out),optional::status
      logical,intent(in),optional::randomize
      real(dp)::flo,fhi,p,u,a,b,c,d
      integer::i,istat
      logical::rnd,disc
      istat=0; rnd=.true.; if(present(randomize))rnd=randomize
      disc=is_discrete_family(result%family)
      allocate(residuals(size(y)))
      if(size(result%mu%fitted)/=size(y))then
         residuals=0.0_dp;istat=1;goto 900
      end if
      do i=1,size(y)
         a=result%mu%fitted(i);b=1.0_dp;c=1.0_dp;d=1.0_dp
         if(allocated(result%sigma%fitted))b=result%sigma%fitted(i)
         if(allocated(result%nu%fitted))c=result%nu%fitted(i)
         if(allocated(result%tau%fitted))d=result%tau%fitted(i)
         if(disc)then
            flo=family_cdf(result%family,y(i)-1.0_dp,a,b,c,d,istat)
            fhi=family_cdf(result%family,y(i),a,b,c,d,istat)
            if(istat/=0)exit
            if(rnd)then;call random_number(u);p=flo+u*(fhi-flo);else;p=0.5_dp*(flo+fhi);end if
         else
            p=family_cdf(result%family,y(i),a,b,c,d,istat)
            if(istat/=0)exit
         end if
         p=min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,p))
         residuals(i)=normal_quantile(p)
      end do
      if(istat/=0)residuals=result%residuals
900   if(present(status))status=istat
   end subroutine randomized_quantile_residuals

   real(dp) function family_cdf(fam,x,a,b,c,d,status) result(p)
      integer,intent(in)::fam
      real(dp),intent(in)::x,a,b,c,d
      integer,intent(out)::status
      status=0
      select case(fam)
      case(GAMLSS_NO); p=pNO(x,a,b)
      case(GAMLSS_GA); p=pGA(x,a,b)
      case(GAMLSS_BE); p=pBE(x,a,b)
      case(GAMLSS_WEI);p=pWEI(x,a,b)
      case(GAMLSS_LNO);p=pLNO(x,a,b,c)
      case(GAMLSS_TF);p=pTF(x,a,b,c)
      case(GAMLSS_EGB2);p=pEGB2(x,a,b,c,d)
      case(GAMLSS_GB2);p=pGB2(x,a,b,c,d)
      case(GAMLSS_GG);p=pGG(x,a,b,c)
      case(GAMLSS_JSUO);p=pJSUo(x,a,b,c,d)
      case(GAMLSS_BCCG);p=pBCCG(x,a,b,c)
      case(GAMLSS_BCT);p=pBCT(x,a,b,c,d)
      case(GAMLSS_BCPE);p=pBCPE(x,a,b,c,d)
      case(GAMLSS_NBI);p=pNBI(x,a,b)
      case(GAMLSS_NBII);p=pNBII(x,a,b)
      case(GAMLSS_ZIP);p=pZIP(x,a,b)
      case(GAMLSS_PIG);p=pPIG(x,a,b)
      case default
         p=0.5_dp;status=2
      end select
   end function family_cdf

   logical pure function is_discrete_family(fam) result(ans)
      integer,intent(in)::fam
      select case(fam)
      case(GAMLSS_NBI,GAMLSS_NBII,GAMLSS_ZIP,GAMLSS_PIG)
         ans=.true.
      case default
         ans=.false.
      end select
   end function is_discrete_family

   subroutine residual_moments(r,mean,sd,skew,excess_kurtosis)
      real(dp),intent(in)::r(:)
      real(dp),intent(out)::mean,sd,skew,excess_kurtosis
      real(dp)::m2,m3,m4,n
      n=real(size(r),dp);mean=sum(r)/max(1.0_dp,n)
      if(size(r)<2)then;sd=0.0_dp;skew=0.0_dp;excess_kurtosis=0.0_dp;return;end if
      m2=sum((r-mean)**2)/n;m3=sum((r-mean)**3)/n;m4=sum((r-mean)**4)/n
      sd=sqrt(max(0.0_dp,m2))
      if(sd>0.0_dp)then;skew=m3/sd**3;excess_kurtosis=m4/sd**4-3.0_dp
      else;skew=0.0_dp;excess_kurtosis=0.0_dp;end if
   end subroutine residual_moments

   subroutine residual_acf(r,max_lag,acf)
      real(dp),intent(in)::r(:)
      integer,intent(in)::max_lag
      real(dp),allocatable,intent(out)::acf(:)
      real(dp)::m,den
      integer::k,n
      n=size(r);allocate(acf(0:max(0,max_lag)));m=sum(r)/real(max(1,n),dp)
      den=sum((r-m)**2);acf(0)=1.0_dp
      do k=1,max_lag
         if(k>=n .or. den<=0.0_dp)then;acf(k)=0.0_dp
         else;acf(k)=sum((r(1:n-k)-m)*(r(1+k:n)-m))/den;end if
      end do
   end subroutine residual_acf

   elemental real(dp) function gaic(global_deviance,df,k) result(v)
      real(dp),intent(in)::global_deviance,df
      real(dp),intent(in),optional::k
      real(dp)::kk
      kk=2.0_dp;if(present(k))kk=k
      v=global_deviance+kk*df
   end function gaic

   elemental real(dp) function likelihood_ratio_stat(reduced_deviance,full_deviance) result(v)
      real(dp),intent(in)::reduced_deviance,full_deviance
      v=max(0.0_dp,reduced_deviance-full_deviance)
   end function likelihood_ratio_stat

   subroutine deviance_increment(reduced,full,increment)
      type(gamlss_result_t),intent(in)::reduced,full
      real(dp),allocatable,intent(out)::increment(:)
      integer::n
      n=min(size(reduced%case_deviance),size(full%case_deviance));allocate(increment(n))
      increment=reduced%case_deviance(1:n)-full%case_deviance(1:n)
   end subroutine deviance_increment

   subroutine hat_values_penalized(x,w,penalty,lambda,hat,status)
      real(dp),intent(in)::x(:,:),w(:),penalty(:,:),lambda
      real(dp),allocatable,intent(out)::hat(:)
      integer,intent(out),optional::status
      real(dp),allocatable::a(:,:),ainv(:,:)
      integer::i,j,k,istat,p
      p=size(x,2);allocate(hat(size(x,1)));hat=0.0_dp
      if(size(w)/=size(x,1).or.size(penalty,1)/=p.or.size(penalty,2)/=p)then
         if(present(status))status=1;return
      end if
      allocate(a(p,p));a=lambda*penalty
      do i=1,size(x,1);do j=1,p;do k=1,p
         a(j,k)=a(j,k)+w(i)*x(i,j)*x(i,k)
      end do;end do;end do
      call invert_matrix(a,ainv,istat)
      if(istat==0)then
         do i=1,size(x,1)
            hat(i)=w(i)*dot_product(x(i,:),matmul(ainv,x(i,:)))
         end do
      end if
      if(present(status))status=istat
   end subroutine hat_values_penalized

   real(dp) function effective_df_penalized(x,w,penalty,lambda) result(edf)
      real(dp),intent(in)::x(:,:),w(:),penalty(:,:),lambda
      real(dp),allocatable::h(:)
      integer::status
      call hat_values_penalized(x,w,penalty,lambda,h,status)
      if(status==0)then;edf=sum(h);else;edf=real(size(x,2),dp);end if
   end function effective_df_penalized

end module gamlss_diagnostics
