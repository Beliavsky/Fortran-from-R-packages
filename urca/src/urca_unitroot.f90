module urca_unitroot
   use urca_kinds, only : dp
   use urca_types, only : ur_test_result, lm_result
   use urca_regression, only : lm_fit, coefficient_t, aic_lm, bic_lm
   implicit none
   private
   integer, parameter, public :: UR_NONE=0, UR_DRIFT=1, UR_TREND=2
   integer, parameter, public :: LAG_FIXED=0, LAG_AIC=1, LAG_BIC=2
   integer, parameter, public :: ERS_DFGLS=1, ERS_PTEST=2
   integer, parameter, public :: PP_ZALPHA=1, PP_ZTAU=2
   integer, parameter, public :: PP_CONSTANT=1, PP_TREND=2
   integer, parameter, public :: KPSS_MU=1, KPSS_TAU=2
   integer, parameter, public :: SP_RHO=1, SP_TAU=2
   integer, parameter, public :: ZA_INTERCEPT=1, ZA_TREND=2, ZA_BOTH=3
   public :: adf_test, ers_test, kpss_test, pp_test, schmidt_phillips_test, zivot_andrews_test
contains
   function diff1(y) result(d)
      real(dp),intent(in)::y(:)
      real(dp),allocatable::d(:)
      allocate(d(max(0,size(y)-1)))
      if(size(y)>1)d=y(2:)-y(:size(y)-1)
   end function diff1

   subroutine adf_design(y,trim,p,model,response,x)
      real(dp),intent(in)::y(:)
      integer,intent(in)::trim,p,model
      real(dp),allocatable,intent(out)::response(:),x(:,:)
      real(dp),allocatable::dy(:)
      integer::nobs,i,j,c
      dy=diff1(y)
      nobs=size(y)-trim-1
      if(nobs<=0)then
      allocate(response(0),x(0,0))
      return
      end if
      select case(model)
      case(UR_NONE)
      c=1+p
      case(UR_DRIFT)
      c=2+p
      case default
      c=3+p
      end select
      allocate(response(nobs),x(nobs,c))
      do i=1,nobs
         response(i)=dy(trim+i)
         select case(model)
         case(UR_NONE)
            x(i,1)=y(trim+i)
            do j=1,p
            x(i,1+j)=dy(trim+i-j)
            end do
         case(UR_DRIFT)
            x(i,1)=1.0_dp
            x(i,2)=y(trim+i)
            do j=1,p
            x(i,2+j)=dy(trim+i-j)
            end do
         case default
            x(i,1)=1.0_dp
            x(i,2)=y(trim+i)
            x(i,3)=real(trim+i,dp)
            do j=1,p
            x(i,3+j)=dy(trim+i-j)
            end do
         end select
      end do
   end subroutine adf_design

   real(dp) function rss_no_predictors(y) result(rss)
      real(dp),intent(in)::y(:)
      rss=sum(y*y)
   end function rss_no_predictors

   function adf_test(y,model,max_lags,selectlags) result(out)
      real(dp),intent(in)::y(:)
      integer,intent(in)::model,max_lags
      integer,intent(in),optional::selectlags
      type(ur_test_result)::out
      type(lm_result)::fit,rf
      real(dp),allocatable::resp(:),x(:,:),xr(:,:),crit(:),cv(:,:)
      real(dp)::bestc,cval,rssr,phi1,phi2,phi3
      integer::sel,p,j,mode,n,rowselec,tauidx,df_r
      if(size(y)<max_lags+4 .or. max_lags<0)then
      out%info=-1
      return
      end if
      mode=LAG_FIXED
      if(present(selectlags))mode=selectlags
      sel=max_lags
      if(mode/=LAG_FIXED .and. max_lags>0)then
         bestc=huge(1.0_dp)
         do p=1,max_lags
            call adf_design(y,max_lags,p,model,resp,x)
            fit=lm_fit(x,resp)
            if(fit%info/=0)cycle
            if(mode==LAG_AIC)cval=aic_lm(fit,2.0_dp)
            if(mode==LAG_BIC)cval=aic_lm(fit,log(real(size(resp),dp)))
            if(cval<bestc)then
            bestc=cval
            sel=p
            end if
         end do
      end if
      call adf_design(y,max_lags,sel,model,resp,x)
      fit=lm_fit(x,resp)
      if(fit%info/=0)then
      out%info=fit%info
      return
      end if
      select case(model)
      case(UR_NONE)
      tauidx=1
      allocate(out%statistic(1))
      case(UR_DRIFT)
      tauidx=2
      allocate(out%statistic(2))
      case default
      tauidx=2
      allocate(out%statistic(3))
      end select
      out%statistic(1)=coefficient_t(fit,tauidx)
      if(model==UR_DRIFT)then
         if(sel>0)then
            allocate(xr(size(resp),sel))
            xr=x(:,3:)
            rf=lm_fit(xr,resp)
            rssr=rf%rss
            df_r=rf%df_resid
         else
         rssr=rss_no_predictors(resp)
         df_r=size(resp)
         end if
         phi1=((rssr-fit%rss)/2.0_dp)/(fit%rss/real(fit%df_resid,dp))
         out%statistic(2)=phi1
      else if(model==UR_TREND)then
         if(sel>0)then
            allocate(xr(size(resp),sel))
            xr=x(:,4:)
            rf=lm_fit(xr,resp)
            rssr=rf%rss
            df_r=rf%df_resid
         else
         rssr=rss_no_predictors(resp)
         df_r=size(resp)
         end if
         phi2=((rssr-fit%rss)/3.0_dp)/(fit%rss/real(fit%df_resid,dp))
         if(allocated(xr))deallocate(xr)
         allocate(xr(size(resp),1+sel))
         xr(:,1)=1.0_dp
         if(sel>0)xr(:,2:)=x(:,4:)
         rf=lm_fit(xr,resp)
         phi3=((rf%rss-fit%rss)/2.0_dp)/(fit%rss/real(fit%df_resid,dp))
         out%statistic(2)=phi2
         out%statistic(3)=phi3
      end if
      n=size(y)-1
      if(n<25)then
      rowselec=1
      else if(n<50)then
      rowselec=2
      else if(n<100)then
      rowselec=3
      else if(n<250)then
      rowselec=4
      else if(n<500)then
      rowselec=5
      else
      rowselec=6
      end if
      select case(model)
      case(UR_NONE)
         allocate(cv(1,3))
         call adf_cv_none(rowselec,cv(1,:))
      case(UR_DRIFT)
         allocate(cv(2,3))
         call adf_cv_drift(rowselec,cv)
      case default
         allocate(cv(3,3))
         call adf_cv_trend(rowselec,cv)
      end select
      out%critical_values=cv
      out%residuals=fit%residuals
      out%coefficients=fit%beta
      allocate(out%std_errors(size(fit%beta)))
      do j=1,size(fit%beta)
      out%std_errors(j)=sqrt(max(0.0_dp,fit%vcov(j,j)))
      end do
      out%lags=sel
      out%info=0
   contains
      subroutine adf_cv_none(r,c)
         integer,intent(in)::r
         real(dp),intent(out)::c(3)
         real(dp),parameter::a(6,3)=reshape([-2.66_dp,-2.62_dp,-2.60_dp,-2.58_dp,-2.58_dp,-2.58_dp,-1.95_dp, &
            & -1.95_dp,-1.95_dp,-1.95_dp,-1.95_dp,-1.95_dp,-1.60_dp,-1.61_dp,-1.61_dp,-1.62_dp,-1.62_dp,-1.62_dp], &
            & [6,3])
         c=a(r,:)
      end subroutine
      subroutine adf_cv_drift(r,cv)
         integer,intent(in)::r
         real(dp),intent(out)::cv(2,3)
         real(dp),parameter::tau(6,3)=reshape([ -3.75_dp,-3.58_dp,-3.51_dp,-3.46_dp,-3.44_dp,-3.43_dp, -3.00_dp, &
            & -2.93_dp,-2.89_dp,-2.88_dp,-2.87_dp,-2.86_dp, -2.63_dp,-2.60_dp,-2.58_dp,-2.57_dp,-2.57_dp,-2.57_dp], &
            & [6,3])
         real(dp),parameter::phi(6,3)=reshape([7.88_dp,7.06_dp,6.70_dp,6.52_dp,6.47_dp,6.43_dp,5.18_dp,4.86_dp, &
            & 4.71_dp,4.63_dp,4.61_dp,4.59_dp,4.12_dp,3.94_dp,3.86_dp,3.81_dp,3.79_dp,3.78_dp],[6,3])
         cv(1,:)=tau(r,:)
         cv(2,:)=phi(r,:)
      end subroutine
      subroutine adf_cv_trend(r,cv)
         integer,intent(in)::r
         real(dp),intent(out)::cv(3,3)
         real(dp),parameter::tau(6,3)=reshape([-4.38_dp,-4.15_dp,-4.04_dp,-3.99_dp,-3.98_dp,-3.96_dp,-3.60_dp, &
            & -3.50_dp,-3.45_dp,-3.43_dp,-3.42_dp,-3.41_dp,-3.24_dp,-3.18_dp,-3.15_dp,-3.13_dp,-3.13_dp,-3.12_dp], &
            & [6,3])
         real(dp),parameter::p2(6,3)=reshape([8.21_dp,7.02_dp,6.50_dp,6.22_dp,6.15_dp,6.09_dp,5.68_dp,5.13_dp, &
            & 4.88_dp,4.75_dp,4.71_dp,4.68_dp,4.67_dp,4.31_dp,4.16_dp,4.07_dp,4.05_dp,4.03_dp],[6,3])
         real(dp),parameter::p3(6,3)=reshape([10.61_dp,9.31_dp,8.73_dp,8.43_dp,8.34_dp,8.27_dp,7.24_dp,6.73_dp, &
            & 6.49_dp,6.49_dp,6.30_dp,6.25_dp,5.91_dp,5.61_dp,5.47_dp,5.47_dp,5.36_dp,5.34_dp],[6,3])
         cv(1,:)=tau(r,:)
         cv(2,:)=p2(r,:)
         cv(3,:)=p3(r,:)
      end subroutine
   end function adf_test

   function kpss_test(y,type,lags_mode,use_lag) result(out)
      real(dp),intent(in)::y(:)
      integer,intent(in)::type
      integer,intent(in),optional::lags_mode,use_lag
      type(ur_test_result)::out
      real(dp),allocatable::res(:),s(:),trend(:),x(:,:)
      type(lm_result)::fit
      real(dp)::nom,s2,den,cov,w
      integer::n,lmax,mode,k
      n=size(y)
      if(n<3)then
      out%info=-1
      return
      end if
      mode=1
      if(present(lags_mode))mode=lags_mode
      if(present(use_lag))then
      lmax=max(0,use_lag)
      else if(mode==2)then
      lmax=int(12.0_dp*(real(n,dp)/100.0_dp)**0.25_dp)
      else if(mode==3)then
      lmax=0
      else
      lmax=int(4.0_dp*(real(n,dp)/100.0_dp)**0.25_dp)
      end if
      allocate(res(n))
      if(type==KPSS_MU)then
         res=y-sum(y)/real(n,dp)
         allocate(out%critical_values(1,4))
         out%critical_values(1,:)=[0.347_dp,0.463_dp,0.574_dp,0.739_dp]
      else
         allocate(x(n,2))
         x(:,1)=1.0_dp
         x(:,2)=[(real(k,dp),k=1,n)]
         fit=lm_fit(x,y)
         res=fit%residuals
         allocate(out%critical_values(1,4))
         out%critical_values(1,:)=[0.119_dp,0.146_dp,0.176_dp,0.216_dp]
      end if
      allocate(s(n))
      s(1)=res(1)
      do k=2,n
      s(k)=s(k-1)+res(k)
      end do
      nom=sum(s*s)/real(n*n,dp)
      s2=sum(res*res)/real(n,dp)
      den=s2
      do k=1,lmax
         cov=dot_product(res(k+1:n),res(1:n-k))
         w=1.0_dp-real(k,dp)/real(lmax+1,dp)
         den=den+2.0_dp*w*cov/real(n,dp)
      end do
      allocate(out%statistic(1))
      out%statistic(1)=nom/den
      out%residuals=res
      out%lags=lmax
   end function kpss_test

   function pp_test(xin,type,model,lags_mode,use_lag) result(out)
      real(dp),intent(in)::xin(:)
      integer,intent(in)::type,model
      integer,intent(in),optional::lags_mode,use_lag
      type(ur_test_result)::out
      integer::n,lmax,mode,k
      real(dp),allocatable::y(:),yl(:),x(:,:),res(:)
      type(lm_result)::fit
      real(dp)::s,myybar,myy,my,mty,sig,lambda,lp,mval,tstat,myt,betat,alpha,cop,w,mean_y,my_stat,beta_stat
      n=size(xin)-1
      if(n<5)then
      out%info=-1
      return
      end if
      allocate(y(n),yl(n))
      y=xin(2:)
      yl=xin(:size(xin)-1)
      mode=1
      if(present(lags_mode))mode=lags_mode
      if(present(use_lag))then
      lmax=max(0,use_lag)
      else if(mode==2)then
      lmax=int(12.0_dp*(real(n,dp)/100.0_dp)**0.25_dp)
      else
      lmax=int(4.0_dp*(real(n,dp)/100.0_dp)**0.25_dp)
      end if
      if(model==PP_TREND)then
         allocate(x(n,3))
         x(:,1)=1
         x(:,2)=yl
         do k=1,n
         x(k,3)=real(k,dp)-real(n,dp)/2.0_dp
         end do
      else
         allocate(x(n,2))
         x(:,1)=1
         x(:,2)=yl
      end if
      fit=lm_fit(x,y)
      if(fit%info/=0)then
      out%info=fit%info
      return
      end if
      res=fit%residuals
      s=sum(res*res)/real(n,dp)
      mean_y=sum(y)/real(n,dp)
      myybar=sum((y-mean_y)**2)/real(n*n,dp)
      myy=sum(y*y)/real(n*n,dp)
      my=sum(y)/(real(n,dp)**1.5_dp)
      sig=s
      do k=1,lmax
      cop=dot_product(res(k+1:n),res(1:n-k))
      w=1.0_dp-real(k,dp)/real(lmax+1,dp)
      sig=sig+2.0_dp*w*cop/real(n,dp)
      end do
      lambda=0.5_dp*(sig-s)
      lp=lambda/sig
      allocate(out%statistic(1),out%critical_values(1,3))
      if(model==PP_TREND)then
         mty=0.0_dp
         do k=1,n
         mty=mty+real(k,dp)*y(k)
         end do
         mty=mty/(real(n,dp)**2.5_dp)
         mval=(1.0_dp-real(n,dp)**(-2))*myy-12*mty*mty+12*(1+1.0_dp/real(n,dp))*mty*my-(4+6.0_dp/real(n, &
            & dp)+2.0_dp/real(n*n,dp))*my*my
         myt=coefficient_t(fit,1)
         betat=coefficient_t(fit,3)
         my_stat=sqrt(s/sig)*myt-lp*sqrt(sig)*my/(sqrt(mval)*sqrt(mval+my*my))
         beta_stat=sqrt(s/sig)*betat-lp*sqrt(sig)*(0.5_dp*my-mty)/(sqrt(mval/12.0_dp)*sqrt(myybar))
         allocate(out%auxiliary_statistics(2))
         out%auxiliary_statistics=[my_stat,beta_stat]
         if(type==PP_ZTAU)then
         tstat=(fit%beta(2)-1.0_dp)/sqrt(fit%vcov(2,2))
         out%statistic(1)=sqrt(s/sig)*tstat-lp*sqrt(sig)/sqrt(mval)
         else
         alpha=fit%beta(2)
         out%statistic(1)=real(n,dp)*(alpha-1.0_dp)-lambda/mval
         end if
         out%critical_values(1,:)=[-3.9638_dp-8.353_dp/n-47.44_dp/(n*n),-3.4126_dp-4.039_dp/n-17.83_dp/(n*n), &
            & -3.1279_dp-2.418_dp/n-7.58_dp/(n*n)]
      else
         myt=coefficient_t(fit,1)
         my_stat=sqrt(s/sig)*myt+lp*sqrt(sig)*my/(sqrt(myy)*sqrt(myybar))
         allocate(out%auxiliary_statistics(1))
         out%auxiliary_statistics(1)=my_stat
         if(type==PP_ZTAU)then
         tstat=(fit%beta(2)-1.0_dp)/sqrt(fit%vcov(2,2))
         out%statistic(1)=sqrt(s/sig)*tstat-lp*sqrt(sig)/sqrt(myybar)
         else
         alpha=fit%beta(2)
         out%statistic(1)=real(n,dp)*(alpha-1.0_dp)-lambda/myybar
         end if
         out%critical_values(1,:)=[-3.4335_dp-5.999_dp/n-29.25_dp/(n*n),-2.8621_dp-2.738_dp/n-8.36_dp/(n*n), &
            & -2.5671_dp-1.438_dp/n-4.48_dp/(n*n)]
      end if
      if(type==PP_ZALPHA)out%critical_values=huge(1.0_dp)
      out%residuals=res
      out%coefficients=fit%beta
      out%lags=lmax
   end function pp_test

   function ers_test(y,type,model,lag_max) result(out)
      real(dp),intent(in)::y(:)
      integer,intent(in)::type,model,lag_max
      type(ur_test_result)::out
      integer::n,k,p,opt,rowsel,i
      real(dp)::ahat,best,bic,sumlc,sig_null,sig_res,what_sq
      real(dp),allocatable::ya(:),z1(:),z2(:),xd(:,:),yd(:),dy(:),resp(:),x(:,:),nullres(:)
      type(lm_result)::fit,fit2
      n=size(y)
      if(n<8 .or. lag_max<0)then
      out%info=-1
      return
      end if
      if(n<50)then
      rowsel=1
      else if(n<100)then
      rowsel=2
      else if(n<=200)then
      rowsel=3
      else
      rowsel=4
      end if
      allocate(ya(n),z1(n),yd(n))
      z1(1)=1.0_dp
      if(model==1)then
         ahat=1.0_dp-7.0_dp/real(n,dp)
         ya(1)=y(1)
         z1(2:)=1.0_dp-ahat
         ya(2:)=y(2:)-ahat*y(:n-1)
         allocate(xd(n,1))
         xd(:,1)=z1
         fit=lm_fit(xd,ya)
         yd=y-fit%beta(1)
      else
         ahat=1.0_dp-13.5_dp/real(n,dp)
         ya(1)=y(1)
         z1(2:)=1.0_dp-ahat
         ya(2:)=y(2:)-ahat*y(:n-1)
         allocate(z2(n),xd(n,2))
         z2(1)=1.0_dp
         do i=2,n
         z2(i)=real(i,dp)-ahat*real(i-1,dp)
         end do
         xd(:,1)=z1
         xd(:,2)=z2
         fit=lm_fit(xd,ya)
         do i=1,n
         yd(i)=y(i)-fit%beta(1)-fit%beta(2)*real(i,dp)
         end do
      end if
      if(type==ERS_DFGLS)then
         call adf_like_no_intercept(yd,lag_max,fit2)
         allocate(out%statistic(1))
         out%statistic(1)=coefficient_t(fit2,1)
         allocate(out%critical_values(1,3))
         if(model==1)then
         out%critical_values(1,:)=[-2.5658_dp-1.960_dp/n-10.04_dp/(n*n),-1.9393_dp-0.398_dp/n,-1.6156_dp-0.181_dp/n]
         else
         call ers_trend_cv(rowsel,out%critical_values(1,:))
         end if
         out%residuals=fit2%residuals
         out%coefficients=fit2%beta
         out%lags=lag_max
      else
         allocate(nullres(n))
         nullres(1)=0
         nullres(2:)=y(2:)-y(:n-1)
         if(model==2)nullres=nullres-sum(nullres)/real(n,dp)
         sig_null=sum(nullres**2)
         sig_res=sum(fit%residuals**2)
         opt=0
         if(lag_max>0)then
            best=huge(1.0_dp)
            do p=1,lag_max
               call adf_level_with_intercept(y,p,fit2)
               bic=bic_lm(fit2)
               if(bic<best)then
               best=bic
               opt=p
               end if
            end do
         else
         call adf_level_with_intercept(y,0,fit2)
         end if
         call adf_level_with_intercept(y,opt,fit2)
         sumlc=0.0_dp
         if(opt>0)sumlc=sum(fit2%beta(3:))
         what_sq=fit2%sigma2/(1.0_dp-sumlc)**2
         allocate(out%statistic(1),out%critical_values(1,3))
         out%statistic(1)=(sig_res-ahat*sig_null)/what_sq
         call ers_p_cv(rowsel,model,out%critical_values(1,:))
         out%residuals=fit%residuals
         out%lags=opt
      end if
      out%detrended=yd
   contains
      subroutine adf_like_no_intercept(z,p,f)
         real(dp),intent(in)::z(:)
         integer,intent(in)::p
         type(lm_result),intent(out)::f
         real(dp),allocatable::dz(:),rr(:),xx(:,:)
         integer::nn,ii,jj
         dz=diff1(z)
         nn=size(z)-p-1
         allocate(rr(nn),xx(nn,1+p))
         do ii=1,nn
         rr(ii)=dz(p+ii)
         xx(ii,1)=z(p+ii)
         do jj=1,p
         xx(ii,1+jj)=dz(p+ii-jj)
         end do
         end do
         f=lm_fit(xx,rr)
      end subroutine
      subroutine adf_level_with_intercept(z,p,f)
         real(dp),intent(in)::z(:)
         integer,intent(in)::p
         type(lm_result),intent(out)::f
         real(dp),allocatable::dz(:),rr(:),xx(:,:)
         integer::nn,ii,jj
         dz=diff1(z)
         nn=size(z)-p-1
         allocate(rr(nn),xx(nn,2+p))
         do ii=1,nn
         rr(ii)=dz(p+ii)
         xx(ii,1)=1
         xx(ii,2)=z(p+ii)
         do jj=1,p
         xx(ii,2+jj)=dz(p+ii-jj)
         end do
         end do
         f=lm_fit(xx,rr)
      end subroutine
      subroutine ers_trend_cv(r,c)
         integer,intent(in)::r
         real(dp),intent(out)::c(3)
         real(dp),parameter::a(4,3)=reshape([-3.77_dp,-3.58_dp,-3.46_dp,-3.48_dp,-3.19_dp,-3.03_dp,-2.93_dp, &
            & -2.89_dp,-2.89_dp,-2.74_dp,-2.64_dp,-2.57_dp],[4,3])
         c=a(r,:)
      end subroutine
      subroutine ers_p_cv(r,m,c)
         integer,intent(in)::r,m
         real(dp),intent(out)::c(3)
         real(dp),parameter::a(4,3,2)=reshape([1.87_dp,1.95_dp,1.91_dp,1.99_dp,2.97_dp,3.11_dp,3.17_dp,3.26_dp, &
            & 3.91_dp,4.17_dp,4.33_dp,4.48_dp,4.22_dp,4.26_dp,4.05_dp,3.96_dp,5.72_dp,5.64_dp,5.66_dp,5.62_dp, &
            & 6.77_dp,6.79_dp,6.86_dp,6.89_dp],[4,3,2])
         c=a(r,:,m)
      end subroutine
   end function ers_test

   function sp_critical(n,type,deg,signif) result(cv)
      integer,intent(in)::n,type,deg
      real(dp),intent(in)::signif
      real(dp)::cv
      real(dp),parameter::obs(7)=[25._dp,50._dp,100._dp,200._dp,500._dp,1000._dp,1e30_dp]
      real(dp),parameter::tau(7,3,4)=reshape(-[3.9_dp,3.73_dp,3.63_dp,3.61_dp,3.59_dp,3.58_dp,3.56_dp,3.18_dp, &
         & 3.11_dp,3.06_dp,3.04_dp,3.04_dp,3.02_dp,3.02_dp,2.85_dp,2.8_dp,2.77_dp,2.76_dp,2.76_dp,2.75_dp,2.75_dp, &
         & 4.52_dp,4.28_dp,4.16_dp,4.12_dp,4.08_dp,4.06_dp,4.06_dp,3.78_dp,3.77_dp,3.65_dp,3.6_dp,3.55_dp,3.55_dp, &
         & 3.53_dp,3.52_dp,3.41_dp,3.34_dp,3.31_dp,3.28_dp,3.26_dp,3.26_dp,3.26_dp,5.07_dp,4.73_dp,4.59_dp,4.53_dp, &
         & 4.5_dp,4.49_dp,4.44_dp,4.26_dp,4.08_dp,4.03_dp,3.99_dp,3.96_dp,3.95_dp,3.93_dp,3.89_dp,3.77_dp,3.72_dp, &
         & 3.69_dp,3.68_dp,3.68_dp,3.67_dp,5.57_dp,5.13_dp,4.99_dp,4.9_dp,4.85_dp,4.83_dp,4.81_dp,4.7_dp,4.47_dp, &
         & 4.39_dp,4.33_dp,4.31_dp,4.31_dp,4.29_dp,4.3_dp,4.15_dp,4.1_dp,4.06_dp,4.03_dp,4.03_dp,4.01_dp],[7,3,4])
      real(dp),parameter::rho(7,3,4)=reshape(-[20.4_dp,22.8_dp,23.8_dp,24.8_dp,25.3_dp,25.3_dp,25.2_dp,15.7_dp, &
         & 17.0_dp,17.5_dp,17.9_dp,18.1_dp,18.1_dp,18.1_dp,13.4_dp,14.3_dp,14.6_dp,14.9_dp,15.0_dp,15.0_dp,15.0_dp, &
         & 24.6_dp,28.4_dp,30.4_dp,31.8_dp,32.4_dp,32.5_dp,32.6_dp,20.1_dp,22.4_dp,23.7_dp,24.2_dp,24.8_dp,24.6_dp, &
         & 24.7_dp,17.8_dp,19.5_dp,20.4_dp,20.7_dp,21.0_dp,21.1_dp,21.1_dp,28.1_dp,33.1_dp,36.3_dp,38.0_dp,39.1_dp, &
         & 39.5_dp,39.7_dp,23.8_dp,27.0_dp,29.1_dp,30.1_dp,30.6_dp,30.8_dp,30.6_dp,21.5_dp,24.0_dp,25.4_dp,26.1_dp, &
         & 26.6_dp,26.7_dp,26.7_dp,31.0_dp,37.4_dp,41.8_dp,44.0_dp,45.3_dp,45.7_dp,45.8_dp,26.9_dp,31.2_dp,34.0_dp, &
         & 35.2_dp,36.2_dp,36.6_dp,36.4_dp,24.7_dp,28.1_dp,30.2_dp,31.2_dp,31.8_dp,32.0_dp,31.9_dp],[7,3,4])
      integer::i,j
      i=7
      do j=1,7
      if(real(n,dp)<=obs(j))then
      i=j
      exit
      end if
      end do
      if(abs(signif-0.05_dp)<1e-12_dp)then
      j=2
      else if(abs(signif-0.1_dp)<1e-12_dp)then
      j=3
      else
      j=1
      end if
      if(type==SP_TAU)then
      cv=tau(i,j,deg)
      else
      cv=rho(i,j,deg)
      end if
   end function sp_critical

   function schmidt_phillips_test(y,type,pol_deg,signif) result(out)
      real(dp),intent(in)::y(:)
      integer,intent(in)::type,pol_deg
      real(dp),intent(in),optional::signif
      type(ur_test_result)::out
      integer::n,lag,i,j,d
      real(dp)::siglev,sigeps,sig,omega2,cop,w,yi,phiy
      real(dp),allocatable::dy(:),trend(:,:),sh(:),xt(:,:),xc(:,:)
      type(lm_result)::fit,cf
      n=size(y)
      d=max(1,min(4,pol_deg))
      if(n<20)then
      out%info=-1
      return
      end if
      siglev=0.01_dp
      if(present(signif))siglev=signif
      lag=int(12.0_dp*(real(n,dp)/100.0_dp)**0.25_dp)
      dy=diff1(y)
      allocate(trend(n,d))
      do j=1,d
      do i=1,n
      trend(i,j)=real(i,dp)**j
      end do
      end do
      allocate(sh(n))
      if(d==1)then
      yi=(y(n)-y(1))/real(n-1,dp)
      phiy=y(1)-yi
      sh=y-phiy-yi*trend(:,1)
      else
         allocate(xt(n-1,d))
         xt(:,1)=1.0_dp
         if(d>2)xt(:,2:d-1)=trend(2:n,1:d-2)
! Regression dy on intercept plus polynomial terms through degree d-1.
         if(d==2)then
         deallocate(xt)
         allocate(xt(n-1,2))
         xt(:,1)=1
         xt(:,2)=trend(2:n,1)
         end if
         if(d==3)then
         deallocate(xt)
         allocate(xt(n-1,3))
         xt(:,1)=1
         xt(:,2:3)=trend(2:n,1:2)
         end if
         if(d==4)then
         deallocate(xt)
         allocate(xt(n-1,4))
         xt(:,1)=1
         xt(:,2:4)=trend(2:n,1:3)
         end if
         fit=lm_fit(xt,dy)
         sh(1)=0
         do i=2,n
         sh(i)=sh(i-1)+fit%residuals(i-1)
         end do
      end if
      if(allocated(xt))deallocate(xt)
      allocate(xt(n-1,d+1))
      xt(:,1)=1
      xt(:,2)=sh(:n-1)
      if(d>1)xt(:,3:)=trend(2:n,1:d-1)
      fit=lm_fit(xt,dy)
      allocate(xc(n-1,d+2))
      xc(:,1)=1
      xc(:,2)=y(:n-1)
      xc(:,3:)=trend(2:n,1:d)
      cf=lm_fit(xc,y(2:))
      sigeps=sum(cf%residuals**2)/real(n,dp)
      sig=sigeps
      do i=1,lag
      if(i>=n-1)exit
      cop=dot_product(cf%residuals(i+1:),cf%residuals(:size(cf%residuals)-i))
      w=(2.0_dp*real(lag-i,dp)/real(lag,dp))**2
      sig=sig+2.0_dp*w*cop/real(n,dp)
      end do
      omega2=sigeps/sig
      allocate(out%statistic(1),out%critical_values(1,1))
      if(type==SP_RHO)then
      out%statistic(1)=real(n,dp)*fit%beta(2)/omega2
      else
      out%statistic(1)=coefficient_t(fit,2)/sqrt(omega2)
      end if
      out%critical_values(1,1)=sp_critical(n,type,d,siglev)
      out%residuals=cf%residuals
      out%coefficients=cf%beta
      out%lags=lag
   end function schmidt_phillips_test

   function zivot_andrews_test(y,model,lag) result(out)
      real(dp),intent(in)::y(:)
      integer,intent(in)::model,lag
      type(ur_test_result)::out
      integer::n,b,i,j,nobs,p,c,bestbp
      real(dp)::stat,best
      real(dp),allocatable::dy(:),resp(:),x(:,:),du(:),dt(:)
      type(lm_result)::fit,bestfit
      n=size(y)
      if(lag<0 .or. n<lag+6)then
      out%info=-1
      return
      end if
      dy=diff1(y)
      nobs=n-lag-1
      allocate(resp(nobs),out%rolling_statistics(n-1))
      out%rolling_statistics=huge(1.0_dp)
      best=huge(1.0_dp)
      bestbp=1
      do b=1,n-1
         select case(model)
         case(ZA_INTERCEPT,ZA_TREND)
         c=4+lag
         case default
         c=5+lag
         end select
         allocate(x(nobs,c),du(n),dt(n))
         du=0
         dt=0
         if(b<n)du(b+1:n)=1
         do i=b+1,n
         dt(i)=real(i-b,dp)
         end do
         do i=1,nobs
            j=lag+i
            resp(i)=y(j+1)
            x(i,1)=1
            x(i,2)=y(j)
            x(i,3)=real(j+1,dp)
            if(lag>0)then
               if(model==ZA_INTERCEPT)then
               do p=1,lag
               x(i,3+p)=dy(j-p)
               end do
               x(i,c)=du(j+1)
               else if(model==ZA_TREND)then
               do p=1,lag
               x(i,3+p)=dy(j-p)
               end do
               x(i,c)=dt(j+1)
               else
               do p=1,lag
               x(i,3+p)=dy(j-p)
               end do
               x(i,c-1)=du(j+1)
               x(i,c)=dt(j+1)
               end if
            else
               if(model==ZA_INTERCEPT)x(i,c)=du(j+1)
               if(model==ZA_TREND)x(i,c)=dt(j+1)
               if(model==ZA_BOTH)then
               x(i,c-1)=du(j+1)
               x(i,c)=dt(j+1)
               end if
            end if
         end do
         fit=lm_fit(x,resp)
         if(fit%info==0 .and. fit%vcov(2,2)>0)then
         stat=(fit%beta(2)-1.0_dp)/sqrt(fit%vcov(2,2))
         out%rolling_statistics(b)=stat
         if(stat<best)then
         best=stat
         bestbp=b
         bestfit=fit
         end if
         end if
         deallocate(x,du,dt)
      end do
      allocate(out%statistic(1),out%critical_values(1,3))
      out%statistic(1)=best
      out%break_point=bestbp
      out%lags=lag
      select case(model)
      case(ZA_INTERCEPT)
      out%critical_values(1,:)=[-5.34_dp,-4.8_dp,-4.58_dp]
      case(ZA_TREND)
      out%critical_values(1,:)=[-4.93_dp,-4.42_dp,-4.11_dp]
      case default
      out%critical_values(1,:)=[-5.57_dp,-5.08_dp,-4.82_dp]
      end select
      if(allocated(bestfit%residuals))then
      out%residuals=bestfit%residuals
      out%coefficients=bestfit%beta
      end if
   end function zivot_andrews_test
end module urca_unitroot
