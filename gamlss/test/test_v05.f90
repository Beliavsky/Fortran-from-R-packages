program test_v05
   use gamlss
   use gamlss_continuous, only : rGA
   use gamlss_continuous_v02, only : rSHASH
   use nlme_types, only : correlation_spec,variance_spec,nlme_control,COR_AR1,VAR_CONSTANT
   implicit none
   call test_correlated_gamma
   call test_multi_random_slopes
   call test_cv_and_full_residuals
   print *, 'test_v05: PASS'
contains
   subroutine assert_true(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then;write(*,'(a)')'FAIL: '//trim(msg);error stop 1;end if
   end subroutine assert_true

   subroutine seed_rng(base)
      integer,intent(in)::base
      integer,allocatable::s(:)
      integer::n,i
      call random_seed(size=n);allocate(s(n));do i=1,n;s(i)=base+53*i;end do;call random_seed(put=s)
   end subroutine seed_rng

   real(dp) function randn() result(z)
      real(dp)::u1,u2
      call random_number(u1);call random_number(u2);u1=max(u1,1.0e-12_dp)
      z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
   end function randn

   subroutine test_correlated_gamma
      integer,parameter::n=96
      real(dp)::y(n),x(n,2),time(n),eta,e,rho
      integer::group(n),i
      type(correlation_spec)::cor
      type(variance_spec)::var
      type(nlme_control)::nctl
      type(gamlss_control_t)::ctl
      type(correlated_rs_result_t)::fit
      call seed_rng(501)
      rho=0.55_dp;e=0.18_dp*randn()
      do i=1,n
         x(i,1)=1.0_dp;x(i,2)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
         time(i)=real(i,dp);group(i)=1
         if(i>1)e=rho*e+0.18_dp*sqrt(1.0_dp-rho*rho)*randn()
         eta=0.45_dp+0.70_dp*x(i,2)+e
         y(i)=rGA(exp(eta),0.30_dp)
      end do
      cor%kind=COR_AR1;cor%fixed=.true.;allocate(cor%par(1));cor%par=rho
      var%kind=VAR_CONSTANT;var%fixed=.true.;allocate(var%par(0))
      nctl=nlme_control();nctl%optimize_covariance=.false.;nctl%reml=.false.
      ctl=gamlss_control_t();ctl%n_cyc=5;ctl%inner_cyc=25;ctl%c_crit=1.0e-4_dp
      call fit_gamlss_correlated_rs(y,x,GAMLSS_GA,fit,correlation=cor,variance=var,time=time,group=group, &
         control=ctl,nlme_control_in=nctl,max_outer=5)
      call assert_true(fit%status==0,'correlated Gamma status')
      call assert_true(abs(fit%model%mu%coefficients(1)-0.45_dp)<0.18_dp,'correlated Gamma intercept')
      call assert_true(abs(fit%model%mu%coefficients(2)-0.70_dp)<0.22_dp,'correlated Gamma slope')
      call assert_true(size(fit%correlation_parameters)==1,'correlated Gamma correlation storage')
      call assert_true(abs(fit%correlation_parameters(1)-rho)<1.0e-12_dp,'fixed AR1 retained')
   end subroutine test_correlated_gamma

   subroutine test_multi_random_slopes
      integer,parameter::ng=8,m=14,n=ng*m
      real(dp)::y(n),xm(n,2),xs(n,1),zr(n,2,2),um(ng,2),us(ng,2),xx,lsig
      integer::group(n),g,i,k
      logical::act(4),corr(4)
      type(multi_random_effects_result_t)::fit
      type(gamlss_control_t)::ctl
      real(dp)::cm,cs
      call seed_rng(702)
      do g=1,ng
         um(g,1)=0.28_dp*sin(0.65_dp*real(g,dp));um(g,2)=0.16_dp*cos(0.51_dp*real(g,dp))
         us(g,1)=0.35_dp*cos(0.58_dp*real(g,dp));us(g,2)=0.20_dp*sin(0.44_dp*real(g,dp))
      end do
      k=0
      do g=1,ng
         do i=1,m
            k=k+1;group(k)=g;xx=-1.0_dp+2.0_dp*real(i-1,dp)/real(m-1,dp)
            xm(k,:)=[1.0_dp,xx];xs(k,1)=1.0_dp
            zr(k,:,1)=[1.0_dp,xx];zr(k,:,2)=[1.0_dp,xx]
            lsig=-1.0_dp+us(g,1)+us(g,2)*xx
            y(k)=1.1_dp+0.62_dp*xx+um(g,1)+um(g,2)*xx+exp(lsig)*randn()
         end do
      end do
      act=[.true.,.true.,.false.,.false.];corr=.true.
      ctl=gamlss_control_t();ctl%n_cyc=18;ctl%inner_cyc=35;ctl%c_crit=2.0e-4_dp
      call fit_gamlss_multi_random_effects(y,xm,zr,group,GAMLSS_NO,fit,active_parameters=act, &
         correlated_within=corr,x_sigma=xs,control=ctl,max_outer=8,tol_cov=2.0e-4_dp)
      call assert_true(fit%status==0,'multi random-slope status')
      call assert_true(all(fit%covariance(1:2,1:2,1)==fit%covariance(1:2,1:2,1)),'mu covariance finite')
      call assert_true(all(fit%covariance(1:2,1:2,2)==fit%covariance(1:2,1:2,2)),'sigma covariance finite')
      call assert_true(fit%covariance(1,1,1)>0.0_dp.and.fit%covariance(2,2,1)>0.0_dp,'mu covariance positive')
      call assert_true(fit%covariance(1,1,2)>0.0_dp.and.fit%covariance(2,2,2)>0.0_dp,'sigma covariance positive')
      cm=correlation(fit%effects(:,1,1),um(:,1));cs=correlation(fit%effects(:,1,2),us(:,1))
      call assert_true(cm>0.45_dp,'multi random mu intercept recovery')
      call assert_true(cs>0.10_dp,'multi random sigma intercept recovery')
   end subroutine test_multi_random_slopes


   subroutine test_cv_and_full_residuals
      integer,parameter::n=72,nr=48
      real(dp)::y(n),xf(n,2),xb(n,1),xx,y2(nr)
      real(dp),allocatable::r(:)
      integer::fold(n),i,status
      type(gamlss_cv_result_t)::cvf,cvb
      type(gamlss_result_t)::fake
      type(gamlss_control_t)::ctl
      call seed_rng(905)
      do i=1,n
         xx=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
         xf(i,:)=[1.0_dp,xx];xb(i,1)=1.0_dp;fold(i)=1+mod(i-1,3)
         y(i)=0.30_dp+1.20_dp*xx+0.38_dp*randn()
      end do
      ctl=gamlss_control_t();ctl%n_cyc=15;ctl%inner_cyc=35;ctl%c_crit=1.0e-4_dp
      call cross_validate_gamlss(y,xf,fold,GAMLSS_NO,cvf,control=ctl)
      call cross_validate_gamlss(y,xb,fold,GAMLSS_NO,cvb,control=ctl)
      call assert_true(cvf%status==0.and.cvb%status==0,'cross-validation status')
      call assert_true(cvf%mean_log_score+0.15_dp<cvb%mean_log_score,'CV prefers true slope model')

      fake%family=GAMLSS_SHASH
      allocate(fake%mu%fitted(nr),fake%sigma%fitted(nr),fake%nu%fitted(nr),fake%tau%fitted(nr))
      fake%mu%fitted=0.2_dp;fake%sigma%fitted=1.1_dp;fake%nu%fitted=0.3_dp;fake%tau%fitted=1.2_dp
      do i=1,nr;y2(i)=rSHASH(0.2_dp,1.1_dp,0.3_dp,1.2_dp);end do
      call randomized_quantile_residuals_all(y2,fake,r,status,randomize=.false.)
      call assert_true(status==0.and.size(r)==nr,'full-family quantile residual status')
      call assert_true(all(r==r).and.maxval(abs(r))<8.0_dp,'full-family quantile residual finite')
   end subroutine test_cv_and_full_residuals

   real(dp) function correlation(a,b) result(c)
      real(dp),intent(in)::a(:),b(:)
      real(dp)::am,bm,da,db
      am=sum(a)/real(size(a),dp);bm=sum(b)/real(size(b),dp);da=sum((a-am)**2);db=sum((b-bm)**2)
      if(da<=0.0_dp.or.db<=0.0_dp)then;c=0.0_dp;else;c=sum((a-am)*(b-bm))/sqrt(da*db);end if
   end function correlation
end program test_v05
