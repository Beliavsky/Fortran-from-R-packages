! Part of the experimental modern Fortran translation of rugarch 1.5-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original rugarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-3.0-only

module rugarch_inference
   use rugarch_kinds,only:dp
   use rugarch_math,only:normal_cdf,regularized_gamma_p
   use rugarch_distributions,only:distribution_cdf
   use rugarch_linalg,only:invert_matrix,covariance_matrix,symmetric_trace_product
   implicit none
   private

   type,public::information_criteria_result
      real(dp)::aic=0.0_dp,bic=0.0_dp,sic=0.0_dp,hqic=0.0_dp
   end type information_criteria_result

   type,public::test_result
      real(dp)::statistic=0.0_dp,p_value=1.0_dp
      real(dp)::parameter1=0.0_dp,parameter2=0.0_dp
      integer::degrees_of_freedom=0
   end type test_result

   type,public::sign_bias_result
      real(dp)::t_value(3)=0.0_dp
      real(dp)::p_value(3)=1.0_dp
      real(dp)::joint_statistic=0.0_dp
      real(dp)::joint_p_value=1.0_dp
      integer::status=0
   end type sign_bias_result

   type,public::nyblom_result
      real(dp),allocatable::individual(:)
      real(dp)::joint=0.0_dp
      real(dp)::critical_10=0.0_dp
      real(dp)::critical_05=0.0_dp
      real(dp)::critical_01=0.0_dp
      integer::status=0
   end type nyblom_result

   type,public::gof_result
      integer,allocatable::groups(:)
      real(dp),allocatable::statistic(:),p_value(:)
   end type gof_result

   abstract interface
      function scalar_objective(x,context) result(value)
         import dp
         real(dp),intent(in)::x(:)
         class(*),intent(in)::context
         real(dp)::value
      end function scalar_objective
   end interface

   public::information_criteria,weighted_box_test,arch_lm_test
   public::sign_bias_test,nyblom_test,adjusted_pearson_gof
   public::numerical_hessian,newey_west_covariance,robust_covariance
   public::classical_covariance

contains

   pure function information_criteria(log_likelihood,nobs,npars) result(ans)
      real(dp),intent(in)::log_likelihood
      integer,intent(in)::nobs,npars
      type(information_criteria_result)::ans
      real(dp)::n,k
      n=real(max(1,nobs),dp);k=real(max(0,npars),dp)
      ans%aic=(-2.0_dp*log_likelihood+2.0_dp*k)/n
      ans%bic=(-2.0_dp*log_likelihood+k*log(n))/n
      ans%sic=(-2.0_dp*log_likelihood)/n+log((n+2.0_dp*k)/n)
      ans%hqic=(-2.0_dp*log_likelihood+2.0_dp*k*log(max(log(n),1.0_dp)))/n
   end function information_criteria

   function weighted_box_test(x,lag,fitdf,weighted,squared,absolute_value,log_squared,monti) result(ans)
      real(dp),intent(in)::x(:)
      integer,intent(in)::lag
      integer,intent(in),optional::fitdf
      logical,intent(in),optional::weighted,squared,absolute_value,log_squared,monti
      type(test_result)::ans
      real(dp),allocatable::y(:),acf(:),pacf(:),obs(:),weights(:)
      real(dp)::mean,den,shape,scale
      integer::n,k,df,j
      logical::use_weighted,use_monti

      n=size(x);k=max(1,min(lag,max(1,n-2)))
      df=0;if(present(fitdf))df=max(0,fitdf)
      use_weighted=.true.;if(present(weighted))use_weighted=weighted
      use_monti=.false.;if(present(monti))use_monti=monti
      allocate(y(n),acf(k),pacf(k),obs(k),weights(k))
      y=x
      if(present(squared))then;if(squared)y=y*y;end if
      if(present(absolute_value))then;if(absolute_value)y=abs(y);end if
      if(present(log_squared))then;if(log_squared)y=log(max(y*y,1.0e-300_dp));end if
      mean=sum(y)/real(n,dp);den=sum((y-mean)**2)
      if(den<=tiny(1.0_dp))return
      do concurrent (j=1:k)
         acf(j)=sum((y(j+1:n)-mean)*(y(1:n-j)-mean))/den
      end do
      call autocorrelation_to_pacf(acf,pacf)
      if(use_monti)then;obs=pacf;else;obs=acf;end if
      if(use_weighted)then
         do concurrent (j=1:k)
            weights(j)=real(k-j+1,dp)/real(k,dp)
         end do
      else
         weights=1.0_dp
      end if
      ans%statistic=real(n,dp)*real(n+2,dp)*sum(weights*obs**2/ &
         real([(n-j,j=1,k)],dp))
      if(use_weighted)then
         shape=0.75_dp*real((k+1)*(k+1)*k,dp)/ &
            real(max(1,2*k*k+3*k+1-6*k*df),dp)
         scale=(2.0_dp/3.0_dp)*real(max(1,2*k*k+3*k+1-6*k*df),dp)/ &
            real(k*(k+1),dp)
         ans%parameter1=shape;ans%parameter2=scale
         ans%p_value=1.0_dp-regularized_gamma_p(shape,ans%statistic/max(scale,1.0e-20_dp))
      else
         ans%degrees_of_freedom=max(1,k-df)
         ans%p_value=chi_square_sf(ans%statistic,real(ans%degrees_of_freedom,dp))
      end if
   end function weighted_box_test

   function arch_lm_test(x,lags,demean) result(ans)
      real(dp),intent(in)::x(:)
      integer,intent(in)::lags
      logical,intent(in),optional::demean
      type(test_result)::ans
      real(dp),allocatable::z(:),design(:,:),target(:),xtx(:,:),inv(:,:),beta(:),fitted(:)
      real(dp)::mean_y,sst,sse
      integer::n,m,j,info
      logical::do_demean

      n=size(x);m=max(1,min(lags,n-3));do_demean=.false.;if(present(demean))do_demean=demean
      allocate(z(n));z=x
      if(do_demean)z=z-sum(z)/real(n,dp)
      z=z*z
      allocate(design(n-m,m+1),target(n-m),xtx(m+1,m+1),inv(m+1,m+1), &
         beta(m+1),fitted(n-m))
      design(:,1)=1.0_dp
      target=z(m+1:n)
      do j=1,m
         design(:,j+1)=z(m+1-j:n-j)
      end do
      xtx=matmul(transpose(design),design)
      call invert_matrix(xtx,inv,info)
      if(info/=0)return
      beta=matmul(inv,matmul(transpose(design),target))
      fitted=matmul(design,beta);mean_y=sum(target)/real(size(target),dp)
      sst=sum((target-mean_y)**2);sse=sum((target-fitted)**2)
      if(sst>0.0_dp)ans%statistic=real(size(target),dp)*max(0.0_dp,1.0_dp-sse/sst)
      ans%degrees_of_freedom=m
      ans%p_value=chi_square_sf(ans%statistic,real(m,dp))
   end function arch_lm_test

   function sign_bias_test(residuals,sigma) result(ans)
      real(dp),intent(in)::residuals(:),sigma(:)
      type(sign_bias_result)::ans
      real(dp),allocatable::z(:),y(:),x(:,:),xtx(:,:),inv(:,:),beta(:),e(:),covb(:,:)
      real(dp)::s2,w(3),joint_matrix(3,3),joint_inverse(3,3)
      integer::n,i,info,df

      n=min(size(residuals),size(sigma));if(n<8)then;ans%status=1;return;end if
      allocate(z(n),y(n-1),x(n-1,4),xtx(4,4),inv(4,4),beta(4),e(n-1),covb(4,4))
      z=residuals(1:n)/max(sigma(1:n),1.0e-20_dp)
      y=z(2:n)**2;x(:,1)=1.0_dp
      do i=1,n-1
         x(i,2)=merge(1.0_dp,0.0_dp,residuals(i)<0.0_dp)
         x(i,3)=x(i,2)*residuals(i)
         x(i,4)=(1.0_dp-x(i,2))*residuals(i)
      end do
      xtx=matmul(transpose(x),x);call invert_matrix(xtx,inv,info)
      if(info/=0)then;ans%status=2;return;end if
      beta=matmul(inv,matmul(transpose(x),y));e=y-matmul(x,beta)
      df=max(1,size(y)-4);s2=sum(e*e)/real(df,dp);covb=s2*inv
      do i=1,3
         ans%t_value(i)=abs(beta(i+1))/sqrt(max(covb(i+1,i+1),1.0e-30_dp))
         ans%p_value(i)=2.0_dp*(1.0_dp-normal_cdf(ans%t_value(i)))
      end do
      w=beta(2:4);joint_matrix=covb(2:4,2:4)
      call invert_matrix(joint_matrix,joint_inverse,info)
      if(info==0)then
         ans%joint_statistic=dot_product(w,matmul(joint_inverse,w))
         ans%joint_p_value=chi_square_sf(ans%joint_statistic,3.0_dp)
      else
         ans%status=3
      end if
   end function sign_bias_test

   function nyblom_test(scores) result(ans)
      real(dp),intent(in)::scores(:,:)
      type(nyblom_result)::ans
      real(dp),allocatable::hes(:,:),hinv(:,:),cum(:,:),xx(:,:)
      integer::n,k,i,j,info

      n=size(scores,1);k=size(scores,2);allocate(ans%individual(k))
      ans%individual=0.0_dp
      if(n<2 .or. k<1)then;ans%status=1;return;end if
      allocate(hes(k,k),hinv(k,k),cum(n,k),xx(k,k))
      hes=matmul(transpose(scores),scores)
      call invert_matrix(hes,hinv,info)
      if(info/=0)then;ans%status=2;return;end if
      cum=0.0_dp
      do j=1,k
         cum(1,j)=scores(1,j)
         do i=2,n
            cum(i,j)=cum(i-1,j)+scores(i,j)
         end do
      end do
      xx=matmul(transpose(cum),cum)
      do j=1,k
         ans%individual(j)=xx(j,j)/max(hes(j,j)*real(n,dp),1.0e-30_dp)
      end do
      ans%joint=symmetric_trace_product(xx,hinv)/real(n,dp)
      call nyblom_critical(k,ans%critical_10,ans%critical_05,ans%critical_01)
   end function nyblom_test

   function adjusted_pearson_gof(z,kind,shape,skew,lambda,groups) result(ans)
      real(dp),intent(in)::z(:),shape,skew,lambda
      integer,intent(in)::kind,groups(:)
      type(gof_result)::ans
      real(dp)::u,expected
      integer::i,j,g,cell
      integer,allocatable::count(:)
      allocate(ans%groups(size(groups)),ans%statistic(size(groups)),ans%p_value(size(groups)))
      ans%groups=groups
      do i=1,size(groups)
         g=max(2,groups(i));allocate(count(g));count=0
         do j=1,size(z)
            u=max(0.0_dp,min(1.0_dp,distribution_cdf(z(j),kind,shape,skew,lambda)))
            cell=min(g,1+int(u*real(g,dp)))
            count(cell)=count(cell)+1
         end do
         expected=real(size(z),dp)/real(g,dp)
         ans%statistic(i)=sum((real(count,dp)-expected)**2/expected)
         ans%p_value(i)=chi_square_sf(ans%statistic(i),real(g-1,dp))
         deallocate(count)
      end do
   end function adjusted_pearson_gof

   subroutine numerical_hessian(fun,x,context,hessian)
      procedure(scalar_objective)::fun
      real(dp),intent(in)::x(:)
      class(*),intent(in)::context
      real(dp),intent(out)::hessian(size(x),size(x))
      real(dp),allocatable::xp(:),xm(:),xpp(:),xpm(:),xmp(:),xmm(:),h(:)
      real(dp)::f0
      integer::n,i,j
      n=size(x);allocate(xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n),h(n))
      h=max(abs(x)*epsilon(1.0_dp)**(1.0_dp/3.0_dp),1.0e-6_dp)
      f0=fun(x,context);hessian=0.0_dp
      do i=1,n
         xp=x;xm=x;xp(i)=xp(i)+h(i);xm(i)=xm(i)-h(i)
         hessian(i,i)=(fun(xp,context)-2.0_dp*f0+fun(xm,context))/(h(i)*h(i))
         do j=i+1,n
            xpp=x;xpm=x;xmp=x;xmm=x
            xpp(i)=xpp(i)+h(i);xpp(j)=xpp(j)+h(j)
            xpm(i)=xpm(i)+h(i);xpm(j)=xpm(j)-h(j)
            xmp(i)=xmp(i)-h(i);xmp(j)=xmp(j)+h(j)
            xmm(i)=xmm(i)-h(i);xmm(j)=xmm(j)-h(j)
            hessian(i,j)=(fun(xpp,context)-fun(xpm,context)-fun(xmp,context)+ &
               fun(xmm,context))/(4.0_dp*h(i)*h(j))
            hessian(j,i)=hessian(i,j)
         end do
      end do
   end subroutine numerical_hessian

   function newey_west_covariance(data,nlag,center) result(cv)
      real(dp),intent(in)::data(:,:)
      integer,intent(in),optional::nlag
      logical,intent(in),optional::center
      real(dp),allocatable::cv(:,:),x(:,:),gamma(:,:)
      integer::n,k,l,bw,j
      logical::do_center
      n=size(data,1);k=size(data,2);allocate(cv(k,k),x(n,k),gamma(k,k));x=data
      do_center=.true.;if(present(center))do_center=center
      if(do_center)then
         do j=1,k;x(:,j)=x(:,j)-sum(x(:,j))/real(max(1,n),dp);end do
      end if
      bw=min(n-1,int(1.2_dp*real(max(1,n),dp)**(1.0_dp/3.0_dp)))
      if(present(nlag))bw=max(0,min(n-1,nlag))
      cv=matmul(transpose(x),x)/real(max(1,n),dp)
      do l=1,bw
         gamma=matmul(transpose(x(l+1:n,:)),x(1:n-l,:))/real(n,dp)
         cv=cv+(1.0_dp-real(l,dp)/real(bw+1,dp))*(gamma+transpose(gamma))
      end do
   end function newey_west_covariance

   subroutine classical_covariance(hessian,nobs,vcv,info)
      real(dp),intent(in)::hessian(:,:)
      integer,intent(in)::nobs
      real(dp),intent(out)::vcv(size(hessian,1),size(hessian,2))
      integer,intent(out)::info
      real(dp),allocatable::a(:,:)
      allocate(a(size(hessian,1),size(hessian,2)))
      a=hessian/real(max(1,nobs),dp)
      call invert_matrix(a,vcv,info)
      if(info==0)vcv=vcv/real(max(1,nobs),dp)
   end subroutine classical_covariance

   subroutine robust_covariance(hessian,scores,nlag,vcv,info)
      real(dp),intent(in)::hessian(:,:),scores(:,:)
      integer,intent(in),optional::nlag
      real(dp),intent(out)::vcv(size(hessian,1),size(hessian,2))
      integer,intent(out)::info
      real(dp),allocatable::a(:,:),ainv(:,:),b(:,:)
      integer::n,lags
      n=size(scores,1);lags=0;if(present(nlag))lags=max(0,nlag)
      allocate(a(size(hessian,1),size(hessian,2)),ainv(size(hessian,1),size(hessian,2)))
      a=hessian/real(max(1,n),dp);call invert_matrix(a,ainv,info)
      if(info/=0)then;vcv=0.0_dp;return;end if
      if(lags>0)then;b=newey_west_covariance(scores,lags,.true.)
      else;b=covariance_matrix(scores,.true.);end if
      vcv=matmul(ainv,matmul(b,ainv))/real(max(1,n),dp)
   end subroutine robust_covariance

   subroutine autocorrelation_to_pacf(acf,pacf)
      real(dp),intent(in)::acf(:)
      real(dp),intent(out)::pacf(size(acf))
      real(dp),allocatable::phi(:,:),v(:)
      integer::k,j
      allocate(phi(size(acf),size(acf)),v(size(acf)));phi=0.0_dp
      if(size(acf)==0)return
      phi(1,1)=acf(1);pacf(1)=phi(1,1);v(1)=max(1.0e-20_dp,1.0_dp-acf(1)**2)
      do k=2,size(acf)
         phi(k,k)=(acf(k)-sum(phi(k-1,1:k-1)*acf(k-1:1:-1)))/v(k-1)
         do j=1,k-1
            phi(k,j)=phi(k-1,j)-phi(k,k)*phi(k-1,k-j)
         end do
         pacf(k)=phi(k,k);v(k)=max(1.0e-20_dp,v(k-1)*(1.0_dp-phi(k,k)**2))
      end do
   end subroutine autocorrelation_to_pacf

   subroutine nyblom_critical(k,c10,c05,c01)
      integer,intent(in)::k
      real(dp),intent(out)::c10,c05,c01
      real(dp),parameter::table(20,3)=reshape([ &
         0.353_dp,0.610_dp,0.846_dp,1.070_dp,1.280_dp,1.490_dp,1.690_dp,1.890_dp,2.100_dp,2.290_dp, &
         2.490_dp,2.690_dp,2.890_dp,3.080_dp,3.260_dp,3.460_dp,3.640_dp,3.830_dp,4.030_dp,4.220_dp, &
         0.470_dp,0.749_dp,1.010_dp,1.240_dp,1.470_dp,1.680_dp,1.900_dp,2.110_dp,2.320_dp,2.540_dp, &
         2.750_dp,2.960_dp,3.150_dp,3.340_dp,3.540_dp,3.750_dp,3.950_dp,4.140_dp,4.330_dp,4.520_dp, &
         0.748_dp,1.070_dp,1.350_dp,1.600_dp,1.880_dp,2.120_dp,2.350_dp,2.590_dp,2.820_dp,3.050_dp, &
         3.270_dp,3.510_dp,3.690_dp,3.900_dp,4.070_dp,4.300_dp,4.510_dp,4.730_dp,4.920_dp,5.130_dp],[20,3])
      integer::i
      i=max(1,min(20,k));c10=table(i,1);c05=table(i,2);c01=table(i,3)
   end subroutine nyblom_critical

   pure elemental function chi_square_sf(x,df) result(value)
      real(dp),intent(in)::x,df
      real(dp)::value
      value=1.0_dp-regularized_gamma_p(0.5_dp*df,0.5_dp*max(x,0.0_dp))
      value=max(0.0_dp,min(1.0_dp,value))
   end function chi_square_sf

end module rugarch_inference
