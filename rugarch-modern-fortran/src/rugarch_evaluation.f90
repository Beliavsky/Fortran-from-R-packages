! Part of the experimental modern Fortran translation of rugarch 1.5-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original rugarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-3.0-only

module rugarch_evaluation
   use rugarch_kinds,only:dp
   use rugarch_math,only:normal_cdf,normal_quantile,regularized_gamma_p
   use rugarch_linalg,only:invert_matrix
   use rugarch_resampling,only:bootstrap_indices,bootstrap_stationary,bootstrap_block
   implicit none
   private

   type,public::var_duration_result
      real(dp)::shape=1.0_dp,scale=1.0_dp
      real(dp)::unrestricted_log_likelihood=0.0_dp
      real(dp)::restricted_log_likelihood=0.0_dp
      real(dp)::lr_statistic=0.0_dp,p_value=1.0_dp
      integer::exceedances=0,status=0
   end type var_duration_result

   type,public::gmm_result
      real(dp)::moment_mean(4)=0.0_dp
      real(dp)::moment_t(4)=0.0_dp
      real(dp)::q_statistic(3)=0.0_dp
      real(dp)::q_p_value(3)=1.0_dp
      real(dp)::joint_statistic=0.0_dp,joint_p_value=1.0_dp
      integer::degrees_of_freedom=0,status=0
   end type gmm_result

   type,public::hong_li_result
      real(dp)::statistic(7)=0.0_dp
      real(dp)::p_value(7)=1.0_dp
   end type hong_li_result

   type,public::mcs_result
      logical,allocatable::included_range(:),included_sq(:)
      real(dp),allocatable::p_value_range(:),p_value_sq(:)
      integer,allocatable::elimination_order_range(:),elimination_order_sq(:)
      integer::status=0
   end type mcs_result

   public::var_duration_test,gmm_test,hong_li_test,mcs_test
   public::loss_squared_error,loss_absolute_error,loss_trading_return
   public::loss_directional_accuracy,loss_asymmetric_absolute
   public::loss_asymmetric_squared,loss_linex,loss_double_linex

contains

   function var_duration_test(alpha,actual,var_forecast) result(ans)
      real(dp),intent(in)::alpha,actual(:),var_forecast(:)
      type(var_duration_result)::ans
      integer,allocatable::hitpos(:),duration(:),censor(:)
      integer::n,nh,i,nd
      real(dp)::lo,hi,b1,b2,f1,f2,best_b,a

      n=min(size(actual),size(var_forecast));nh=count(actual(1:n)<var_forecast(1:n))
      ans%exceedances=nh
      if(nh<1)then;ans%status=1;return;end if
      allocate(hitpos(nh));hitpos=pack([(i,i=1,n)],actual(1:n)<var_forecast(1:n))
      nd=max(1,nh-1)+merge(1,0,hitpos(1)>1)+merge(1,0,hitpos(nh)<n)
      allocate(duration(nd),censor(nd));duration=0;censor=0;i=0
      if(hitpos(1)>1)then;i=i+1;duration(i)=hitpos(1);censor(i)=1;end if
      if(nh>1)then
         duration(i+1:i+nh-1)=hitpos(2:nh)-hitpos(1:nh-1)
         i=i+nh-1
      end if
      if(hitpos(nh)<n)then;i=i+1;duration(i)=n-hitpos(nh);censor(i)=1;end if
      if(i==0)then;i=1;duration(1)=n;censor(1)=1;end if

      lo=0.001_dp;hi=10.0_dp
      do i=1,100
         b1=lo+(hi-lo)/3.0_dp;b2=hi-(hi-lo)/3.0_dp
         f1=duration_negloglik(b1,duration(1:nd),censor(1:nd))
         f2=duration_negloglik(b2,duration(1:nd),censor(1:nd))
         if(f1<f2)then;hi=b2;else;lo=b1;end if
      end do
      best_b=0.5_dp*(lo+hi)
      ans%shape=best_b
      a=duration_scale(best_b,duration(1:nd),censor(1:nd));ans%scale=1.0_dp/max(a,1.0e-30_dp)
      ans%unrestricted_log_likelihood=-duration_negloglik(best_b,duration(1:nd),censor(1:nd))
      ans%restricted_log_likelihood=-duration_negloglik(1.0_dp,duration(1:nd),censor(1:nd))
      ans%lr_statistic=max(0.0_dp,2.0_dp*(ans%unrestricted_log_likelihood-ans%restricted_log_likelihood))
      ans%p_value=chi_square_sf(ans%lr_statistic,1.0_dp)
      if(alpha<=0.0_dp .or. alpha>=1.0_dp)ans%status=2
   end function var_duration_test

   function gmm_test(z,lags,skew,kurtosis) result(ans)
      real(dp),intent(in)::z(:)
      integer,intent(in),optional::lags
      real(dp),intent(in),optional::skew,kurtosis
      type(gmm_result)::ans
      real(dp),allocatable::m(:,:),g(:),s(:,:),sinv(:,:),block(:,:),gb(:),sb(:,:),sbinv(:,:)
      real(dp)::sk,ku,variance
      integer::n,l,nobs,k,i,j,idx,info

      l=1;if(present(lags))l=max(1,lags);sk=0.0_dp;if(present(skew))sk=skew
      ku=3.0_dp;if(present(kurtosis))ku=kurtosis
      n=size(z);nobs=n-l;k=4+3*l
      if(nobs<=k+1)then;ans%status=1;return;end if
      allocate(m(nobs,k),g(k),s(k,k),sinv(k,k));m=0.0_dp
      do i=1,nobs
         m(i,1)=z(i+l)
         m(i,2)=z(i+l)**2-1.0_dp
         m(i,3)=z(i+l)**3-sk
         m(i,4)=z(i+l)**4-ku
         idx=5
         do j=1,l;m(i,idx)=m(i,2)*(z(i+l-j)**2-1.0_dp);idx=idx+1;end do
         do j=1,l;m(i,idx)=m(i,3)*(z(i+l-j)**3-sk);idx=idx+1;end do
         do j=1,l;m(i,idx)=m(i,4)*(z(i+l-j)**4-ku);idx=idx+1;end do
      end do
      g=sum(m,dim=1)/real(nobs,dp)
      do i=1,4
         ans%moment_mean(i)=g(i);variance=sum((m(:,i)-g(i))**2)/real(nobs-1,dp)
         ans%moment_t(i)=g(i)/sqrt(max(variance/real(nobs,dp),1.0e-30_dp))
      end do
      s=matmul(transpose(m),m)/real(nobs,dp);call invert_matrix(s,sinv,info)
      if(info/=0)then;ans%status=2;return;end if
      ans%joint_statistic=real(nobs,dp)*dot_product(g,matmul(sinv,g))
      ans%degrees_of_freedom=k;ans%joint_p_value=chi_square_sf(ans%joint_statistic,real(k,dp))
      allocate(block(nobs,l),gb(l),sb(l,l),sbinv(l,l))
      do idx=1,3
         block=m(:,5+(idx-1)*l:4+idx*l);gb=sum(block,dim=1)/real(nobs,dp)
         sb=matmul(transpose(block),block)/real(nobs,dp);call invert_matrix(sb,sbinv,info)
         if(info==0)then
            ans%q_statistic(idx)=real(nobs,dp)*dot_product(gb,matmul(sbinv,gb))
            ans%q_p_value(idx)=chi_square_sf(ans%q_statistic(idx),real(l,dp))
         else
            ans%status=max(ans%status,3)
         end if
      end do
   end function gmm_test

   function hong_li_test(pit,lags) result(ans)
      real(dp),intent(in)::pit(:)
      integer,intent(in),optional::lags
      type(hong_li_result)::ans
      integer::p,t,i
      real(dp)::hpit,acon2,acon11,aconm1,vcon,qsum
      p=4;if(present(lags))p=max(1,lags);t=size(pit)
      if(t<20)then;ans%p_value=1.0_dp;return;end if
      ans%statistic(1)=moment_stat(1,1,p,pit)
      ans%statistic(2)=moment_stat(2,2,p,pit)
      ans%statistic(3)=moment_stat(3,3,p,pit)
      ans%statistic(4)=moment_stat(4,4,p,pit)
      ans%statistic(5)=moment_stat(1,2,p,pit)
      ans%statistic(6)=moment_stat(2,1,p,pit)
      hpit=sample_sd(pit)*real(t,dp)**(-1.0_dp/6.0_dp)
      hpit=max(0.01_dp,min(0.25_dp,hpit))
      acon2=simpson_funb(400)
      acon11=(1.0_dp/hpit-2.0_dp)*(5.0_dp/7.0_dp)
      aconm1=(acon11+2.0_dp*acon2)**2-1.0_dp
      vcon=2.0_dp*(50.0_dp/49.0_dp-300.0_dp/294.0_dp+1950.0_dp/1960.0_dp- &
         900.0_dp/1568.0_dp+450.0_dp/2304.0_dp)**2
      qsum=0.0_dp
      do i=1,p
         qsum=qsum+(real(t-i,dp)*hpit*ghat_integral(pit,i,hpit)-hpit*aconm1)/sqrt(vcon)
      end do
      ans%statistic(7)=qsum/sqrt(real(p,dp))
      do i=1,7
         ans%p_value(i)=1.0_dp-normal_cdf(ans%statistic(i))
      end do
   end function hong_li_test

   function mcs_test(losses,alpha,nboot,nblock,method) result(ans)
      real(dp),intent(in)::losses(:,:),alpha
      integer,intent(in),optional::nboot,nblock,method
      type(mcs_result)::ans
      integer::nb,bl,boot_method,n,m,b,i,j,r,ninc,remove,idx
      integer,allocatable::bs(:,:),order(:)
      logical,allocatable::included(:)
      real(dp),allocatable::dbar(:,:),dstar(:,:,:),variance(:,:),z0(:,:,:),zdata(:,:), &
         pvals(:),dibar(:),vardi(:),tx(:),emp(:)
      real(dp)::tr,tsq,maxp

      n=size(losses,1);m=size(losses,2);nb=100;if(present(nboot))nb=max(20,nboot)
      bl=1;if(present(nblock))bl=max(1,nblock);boot_method=bootstrap_stationary
      if(present(method))boot_method=method
      allocate(ans%included_range(m),ans%included_sq(m),ans%p_value_range(m),ans%p_value_sq(m), &
         ans%elimination_order_range(m),ans%elimination_order_sq(m),bs(n,nb))
      call bootstrap_indices(n,nb,bl,boot_method,bs)
      allocate(dbar(m,m),dstar(m,m,nb),variance(m,m),z0(m,m,nb),zdata(m,m))
      do j=1,m
         do i=1,m
            dbar(j,i)=sum(losses(:,i)-losses(:,j))/real(n,dp)
         end do
      end do
      do b=1,nb
         do j=1,m
            do i=1,m
               dstar(j,i,b)=sum(losses(bs(:,b),i)-losses(bs(:,b),j))/real(n,dp)
            end do
         end do
      end do
      variance=0.0_dp
      do b=1,nb;variance=variance+(dstar(:,:,b)-dbar)**2;end do
      variance=variance/real(nb,dp)
      do i=1,m;variance(i,i)=variance(i,i)+1.0_dp;end do
      do b=1,nb;z0(:,:,b)=(dstar(:,:,b)-dbar)/sqrt(max(variance,1.0e-30_dp));end do
      zdata=dbar/sqrt(max(variance,1.0e-30_dp))
      allocate(pvals(m),order(m),included(m),dibar(m),vardi(m),tx(m),emp(nb))

      call mcs_eliminate(.true.)
      ans%p_value_range=pvals;ans%elimination_order_range=order
      idx=1;do while(idx<=m .and. pvals(idx)<alpha);idx=idx+1;end do
      ans%included_range=.false.;if(idx<=m)ans%included_range(order(idx:m))=.true.

      call mcs_eliminate(.false.)
      ans%p_value_sq=pvals;ans%elimination_order_sq=order
      idx=1;do while(idx<=m .and. pvals(idx)<alpha);idx=idx+1;end do
      ans%included_sq=.false.;if(idx<=m)ans%included_sq(order(idx:m))=.true.

   contains
      subroutine mcs_eliminate(use_range)
         logical,intent(in)::use_range
         included=.true.;pvals=1.0_dp;order=0
         do r=1,m-1
            ninc=count(included)
            if(use_range)then
               tr=maxval(zdata,mask=spread(included,1,m).and.spread(included,2,m))
               do b=1,nb
                  emp(b)=maxval(abs(z0(:,:,b)),mask=spread(included,1,m).and.spread(included,2,m))
               end do
            else
               tsq=0.5_dp*sum(zdata*zdata,mask=spread(included,1,m).and.spread(included,2,m))
               tr=tsq
               do b=1,nb
                  emp(b)=0.5_dp*sum(z0(:,:,b)**2,mask=spread(included,1,m).and.spread(included,2,m))
               end do
            end if
            pvals(r)=real(count(emp>tr),dp)/real(nb,dp)
            dibar=0.0_dp;vardi=0.0_dp;tx=-huge(1.0_dp)
            do i=1,m
               if(.not.included(i))cycle
               do j=1,m
                  if(included(j))dibar(i)=dibar(i)+dbar(j,i)
               end do
               dibar(i)=dibar(i)/real(max(1,ninc-1),dp)
               do b=1,nb
                  tr=0.0_dp
                  do j=1,m
                     if(included(j))tr=tr+dstar(j,i,b)
                  end do
                  tr=tr/real(max(1,ninc-1),dp)
                  vardi(i)=vardi(i)+(tr-dibar(i))**2
               end do
               vardi(i)=vardi(i)/real(nb,dp)
               tx(i)=dibar(i)/sqrt(max(vardi(i),1.0e-30_dp))
            end do
            remove=maxloc(tx,dim=1);order(r)=remove;included(remove)=.false.
         end do
         order(m)=maxloc(merge(1.0_dp,0.0_dp,included),dim=1)
         maxp=pvals(1)
         do i=2,m
            maxp=max(maxp,pvals(i));pvals(i)=maxp
         end do
      end subroutine mcs_eliminate
   end function mcs_test

   pure elemental function loss_squared_error(actual,forecast) result(value)
      real(dp),intent(in)::actual,forecast;real(dp)::value;value=(actual-forecast)**2
   end function loss_squared_error
   pure elemental function loss_absolute_error(actual,forecast) result(value)
      real(dp),intent(in)::actual,forecast;real(dp)::value;value=abs(actual-forecast)
   end function loss_absolute_error
   pure elemental function loss_trading_return(actual,forecast) result(value)
      real(dp),intent(in)::actual,forecast;real(dp)::value;value=-sign(1.0_dp,forecast)*actual
   end function loss_trading_return
   pure elemental function loss_directional_accuracy(actual,forecast) result(value)
      real(dp),intent(in)::actual,forecast;real(dp)::value
      value=-sign(1.0_dp,forecast)*sign(1.0_dp,actual)
   end function loss_directional_accuracy
   pure elemental function loss_asymmetric_absolute(actual,forecast,alpha) result(value)
      real(dp),intent(in)::actual,forecast,alpha;real(dp)::value,e
      e=actual-forecast;value=(alpha-merge(1.0_dp,0.0_dp,actual<forecast))*e
   end function loss_asymmetric_absolute
   pure elemental function loss_asymmetric_squared(actual,forecast,alpha) result(value)
      real(dp),intent(in)::actual,forecast,alpha;real(dp)::value,e
      e=actual-forecast;value=abs(alpha-merge(1.0_dp,0.0_dp,actual<forecast))*e*e
   end function loss_asymmetric_squared
   pure elemental function loss_linex(actual,forecast,alpha) result(value)
      real(dp),intent(in)::actual,forecast,alpha;real(dp)::value,e
      e=actual-forecast;value=exp(alpha*e)-alpha*e-1.0_dp
   end function loss_linex
   pure elemental function loss_double_linex(actual,forecast,alpha,beta) result(value)
      real(dp),intent(in)::actual,forecast,alpha,beta;real(dp)::value,e
      e=actual-forecast;value=exp(alpha*e)+exp(-beta*e)-(alpha-beta)*e-2.0_dp
   end function loss_double_linex

   pure function duration_scale(b,d,c) result(a)
      real(dp),intent(in)::b;integer,intent(in)::d(:),c(:);real(dp)::a
      integer::events
      events=size(d)-c(1)-c(size(c))
      a=(real(max(1,events),dp)/sum(real(d,dp)**b))**(1.0_dp/b)
   end function duration_scale

   pure function duration_negloglik(b,d,c) result(value)
      real(dp),intent(in)::b;integer,intent(in)::d(:),c(:);real(dp)::value,a,ld
      integer::i
      if(b<=0.0_dp)then;value=huge(1.0_dp);return;end if
      a=duration_scale(b,d,c);value=0.0_dp
      do i=1,size(d)
         if(c(i)==1)then
            ld=-(a*real(d(i),dp))**b
         else
            ld=b*log(a)+log(b)+(b-1.0_dp)*log(real(d(i),dp))-(a*real(d(i),dp))**b
         end if
         value=value-ld
      end do
   end function duration_negloglik

   pure function moment_stat(m,l,p,pit) result(value)
      integer,intent(in)::m,l,p;real(dp),intent(in)::pit(:);real(dp)::value
      real(dp)::part1,part2,part3,w,r0,rj
      integer::j
      part1=0.0_dp;part2=0.0_dp;part3=0.0_dp;r0=rcc(0,m,l,pit)
      do j=1,p
         w=1.0_dp-real(j,dp)/real(p,dp);rj=rcc(j,m,l,pit)
         part1=part1+w*w*real(size(pit)-j,dp)*(rj/max(abs(r0),1.0e-20_dp))**2
         part2=part2+w*w;part3=part3+w**4
      end do
      value=(part1-part2)/max(part3,1.0e-20_dp)
   end function moment_stat

   pure function rcc(j,m,l,pit) result(value)
      integer,intent(in)::j,m,l;real(dp),intent(in)::pit(:);real(dp)::value,r1,r2,r3
      integer::t
      t=size(pit)
      if(j==0)then
         r1=sum(pit**(m+l))/real(t,dp);r2=sum(pit**m)/real(t,dp);r3=sum(pit**l)/real(t,dp)
      else
         r1=sum(pit(j+1:t)**m*pit(1:t-j)**l)/real(t,dp)
         r2=sum(pit(j+1:t)**m)/real(t,dp);r3=sum(pit(1:t-j)**l)/real(t,dp)
      end if
      value=r1-r2*r3
   end function rcc

   function ghat_integral(pit,lag,h) result(value)
      real(dp),intent(in)::pit(:),h;integer,intent(in)::lag;real(dp)::value
      real(dp),parameter::nodes(12)=[-0.9815606_dp,-0.9041173_dp,-0.7699027_dp,-0.5873180_dp, &
         -0.3678315_dp,-0.1252334_dp,0.1252334_dp,0.3678315_dp,0.5873180_dp,0.7699027_dp, &
         0.9041173_dp,0.9815606_dp]
      real(dp),parameter::weights(12)=[0.04717534_dp,0.10693933_dp,0.16007833_dp,0.20316743_dp, &
         0.23349254_dp,0.24914705_dp,0.24914705_dp,0.23349254_dp,0.20316743_dp,0.16007833_dp, &
         0.10693933_dp,0.04717534_dp]
      real(dp)::x,y,g
      integer::i,j
      value=0.0_dp
      do i=1,12
         x=0.5_dp*(nodes(i)+1.0_dp)
         do j=1,12
            y=0.5_dp*(nodes(j)+1.0_dp);g=kernel_density_product(x,y,pit,lag,h)
            value=value+0.25_dp*weights(i)*weights(j)*(g-1.0_dp)**2
         end do
      end do
   end function ghat_integral

   function kernel_density_product(x,y,pit,lag,h) result(value)
      real(dp),intent(in)::x,y,pit(:),h;integer,intent(in)::lag;real(dp)::value
      real(dp)::b1,b2
      integer::t,n
      n=size(pit);value=0.0_dp
      do t=1,n-lag
         b1=bounded_kernel(x,pit(t+lag),h);b2=bounded_kernel(y,pit(t),h)
         value=value+b1*b2
      end do
      value=value/real(n-lag,dp)
   end function kernel_density_product

   function bounded_kernel(x,y,h) result(value)
      real(dp),intent(in)::x,y,h;real(dp)::value,den,lo,hi
      value=quartic((x-y)/h)/h
      if(x<h)then
         lo=-x/h;hi=1.0_dp;den=quartic_integral(lo,hi);value=value/max(den,1.0e-20_dp)
      else if(x>1.0_dp-h)then
         lo=-1.0_dp;hi=(1.0_dp-x)/h;den=quartic_integral(lo,hi);value=value/max(den,1.0e-20_dp)
      end if
   end function bounded_kernel

   pure elemental function quartic(z) result(value)
      real(dp),intent(in)::z;real(dp)::value
      if(abs(z)<=1.0_dp)then;value=(15.0_dp/16.0_dp)*(1.0_dp-z*z)**2;else;value=0.0_dp;end if
   end function quartic

   pure function quartic_integral(a,b) result(value)
      real(dp),intent(in)::a,b;real(dp)::value
      value=(15.0_dp/16.0_dp)*((b-a)-(2.0_dp/3.0_dp)*(b**3-a**3)+(1.0_dp/5.0_dp)*(b**5-a**5))
   end function quartic_integral

   function simpson_funb(nstep) result(value)
      integer,intent(in)::nstep;real(dp)::value,h,x
      integer::i,n
      n=nstep;if(mod(n,2)/=0)n=n+1;h=1.0_dp/real(n,dp);value=funb(0.0_dp)+funb(1.0_dp)
      do i=1,n-1;x=real(i,dp)*h;value=value+merge(4.0_dp,2.0_dp,mod(i,2)==1)*funb(x);end do
      value=value*h/3.0_dp
   end function simpson_funb

   pure function funb(b) result(value)
      real(dp),intent(in)::b;real(dp)::value,tmp1,tmp2
      tmp1=(8.0_dp/15.0_dp+b-(2.0_dp/3.0_dp)*b**3+(1.0_dp/5.0_dp)*b**5)**(-2)
      tmp2=b*(1.0_dp-b*b)**4+128.0_dp/315.0_dp+(8.0_dp/3.0_dp)*b**3- &
         (24.0_dp/5.0_dp)*b**5+(24.0_dp/7.0_dp)*b**7-(8.0_dp/9.0_dp)*b**9
      value=tmp1*tmp2
   end function funb

   pure function sample_sd(x) result(value)
      real(dp),intent(in)::x(:);real(dp)::value,m
      m=sum(x)/real(max(1,size(x)),dp)
      value=sqrt(sum((x-m)**2)/real(max(1,size(x)-1),dp))
   end function sample_sd

   pure elemental function chi_square_sf(x,df) result(value)
      real(dp),intent(in)::x,df;real(dp)::value
      value=1.0_dp-regularized_gamma_p(0.5_dp*df,0.5_dp*max(x,0.0_dp))
      value=max(0.0_dp,min(1.0_dp,value))
   end function chi_square_sf

end module rugarch_evaluation
