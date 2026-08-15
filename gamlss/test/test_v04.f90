program test_v04
   use gamlss
   use nlme_types, only : correlation_spec,variance_spec,nlme_control,COR_AR1,VAR_CONSTANT
   use gamlss_special, only : normal_quantile
   implicit none
   call test_correlated_no
   call test_multi_random
   call test_sigma_selection
   call test_diagnostics
   print *, 'test_v04: PASS'
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
      call random_seed(size=n);allocate(s(n))
      do i=1,n;s(i)=base+37*i;end do
      call random_seed(put=s)
   end subroutine seed_rng

   real(dp) function randn() result(z)
      real(dp)::u1,u2
      call random_number(u1);call random_number(u2)
      u1=max(u1,1.0e-12_dp);z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
   end function randn

   subroutine test_correlated_no
      integer,parameter::n=140
      real(dp)::y(n),x(n,2),time(n),e,rho,sd
      integer::group(n),i
      type(correlation_spec)::cor
      type(variance_spec)::var
      type(nlme_control)::ctl
      type(correlated_no_result_t)::fit
      call seed_rng(104)
      rho=0.55_dp;sd=0.45_dp;e=sd*randn()
      do i=1,n
         x(i,1)=1.0_dp;x(i,2)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
         time(i)=real(i,dp);group(i)=1
         if(i>1)e=rho*e+sd*sqrt(1.0_dp-rho*rho)*randn()
         y(i)=0.8_dp+1.1_dp*x(i,2)+e
      end do
      cor%kind=COR_AR1;allocate(cor%par(1));cor%par=0.2_dp
      var%kind=VAR_CONSTANT;allocate(var%par(0))
      ctl=nlme_control();ctl%reml=.false.;ctl%max_iter=160;ctl%tolerance=1.0e-6_dp
      call fit_gamlss_no_gls(y,x,fit,correlation=cor,variance=var,time=time,group=group,control=ctl)
      call assert_true(fit%status==0,'correlated NO fit status')
      call assert_true(abs(fit%gls%beta(1)-0.8_dp)<0.18_dp,'correlated NO intercept')
      call assert_true(abs(fit%gls%beta(2)-1.1_dp)<0.18_dp,'correlated NO slope')
      call assert_true(allocated(fit%gls%correlation_parameters),'AR1 parameter allocation')
      call assert_true(abs(fit%gls%correlation_parameters(1)-rho)<0.22_dp,'AR1 recovery')
   end subroutine test_correlated_no

   subroutine test_multi_random
      integer,parameter::ng=10,m=20,n=ng*m
      real(dp)::y(n),xm(n,2),xs(n,1),u_mu(ng),u_sig(ng),xx,sd
      integer::g(n),i,j,k
      logical::active(4)
      type(multi_random_intercept_result_t)::fit
      type(gamlss_control_t)::ctl
      real(dp)::cm,cs
      call seed_rng(208)
      do j=1,ng
         u_mu(j)=0.30_dp*sin(0.8_dp*real(j,dp))
         u_sig(j)=0.22_dp*cos(0.7_dp*real(j,dp))
      end do
      k=0
      do j=1,ng
         do i=1,m
            k=k+1;g(k)=j;xx=-1.0_dp+2.0_dp*real(i-1,dp)/real(m-1,dp)
            xm(k,:)=[1.0_dp,xx];xs(k,1)=1.0_dp;sd=exp(-0.9_dp+u_sig(j))
            y(k)=1.2_dp+0.65_dp*xx+u_mu(j)+sd*randn()
         end do
      end do
      active=[.true.,.true.,.false.,.false.]
      ctl=gamlss_control_t();ctl%n_cyc=35;ctl%inner_cyc=70;ctl%c_crit=2.0e-4_dp
      call fit_gamlss_multi_random_intercept(y,xm,g,GAMLSS_NO,fit,active_parameters=active,x_sigma=xs,control=ctl)
      call assert_true(fit%status==0,'multi-parameter random-intercept status')
      call assert_true(fit%lambda(1)>0.0_dp.and.fit%lambda(2)>0.0_dp,'multi random lambdas')
      cm=correlation(fit%effects(:,1),u_mu);cs=correlation(fit%effects(:,2),u_sig)
      call assert_true(cm>0.55_dp,'mu random-effect recovery')
      call assert_true(cs>0.35_dp,'sigma random-effect recovery')
   end subroutine test_multi_random

   subroutine test_sigma_selection
      integer,parameter::n=180
      real(dp)::y(n),xmu(n,1),base(n,1),cand(n,2),z1,z2,sig
      integer::i
      type(stepwise_result_t)::sel
      type(gamlss_control_t)::ctl
      call seed_rng(311)
      do i=1,n
         z1=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
         z2=sin(0.37_dp*real(i,dp))
         xmu(i,1)=1.0_dp;base(i,1)=1.0_dp;cand(i,:)=[z1,z2]
         sig=exp(-0.65_dp+0.75_dp*z1);y(i)=0.4_dp+sig*randn()
      end do
      ctl=gamlss_control_t();ctl%n_cyc=30
      call stepwise_gaic_parameter(y,base,cand,GAMLSS_NO,2,sel,direction=STEP_BOTH, &
         k=log(real(n,dp)),x_mu=xmu,control=ctl)
      call assert_true(sel%status==0,'sigma stepwise status')
      call assert_true(size(sel%selected)==1,'sigma stepwise selected count')
      call assert_true(sel%selected(1)==1,'sigma stepwise selected true predictor')
   end subroutine test_sigma_selection

   subroutine test_diagnostics
      integer,parameter::n=120
      real(dp)::r(n),x(n),jb
      real(dp),allocatable::hat(:),stud(:),cook(:)
      type(worm_result_t)::wr
      integer::i,status
      do i=1,n
         x(i)=real(mod(37*i,n),dp);r(i)=normal_quantile((real(i,dp)-0.375_dp)/(real(n,dp)+0.25_dp))
      end do
      call worm_plot_diagnostics(r,x,3,wr,status)
      call assert_true(status==0,'worm diagnostic status')
      call assert_true(maxval(abs(wr%coefficients))<0.30_dp,'worm coefficients near zero')
      allocate(hat(n));hat=0.05_dp
      call influence_from_hat(r,hat,2,stud,cook,status)
      call assert_true(status==0.and.all(cook>=0.0_dp),'influence diagnostics')
      jb=jarque_bera_statistic(r)
      call assert_true(jb<2.0_dp,'normal-score Jarque-Bera')
   end subroutine test_diagnostics

   real(dp) function correlation(a,b) result(c)
      real(dp),intent(in)::a(:),b(:)
      real(dp)::am,bm,da,db
      am=sum(a)/real(size(a),dp);bm=sum(b)/real(size(b),dp)
      da=sum((a-am)**2);db=sum((b-bm)**2)
      if(da<=0.0_dp.or.db<=0.0_dp)then;c=0.0_dp;else;c=sum((a-am)*(b-bm))/sqrt(da*db);end if
   end function correlation
end program test_v04
