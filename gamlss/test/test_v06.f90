program test_v06
   use gamlss
   use gamlss_continuous, only : qGA
   use gamlss_special, only : normal_cdf
   use nlme_types, only : correlation_spec,COR_AR1
   implicit none
   call test_gaussian_copula_gamma
   call test_cross_parameter_random_covariance
   print *, 'test_v06: PASS'
contains
   subroutine assert_true(ok,msg)
      logical,intent(in) :: ok
      character(*),intent(in) :: msg
      if(.not.ok)then;write(*,'(a)')'FAIL: '//trim(msg);error stop 1;end if
   end subroutine assert_true

   subroutine seed_rng(base)
      integer,intent(in) :: base
      integer,allocatable :: s(:)
      integer :: n,i
      call random_seed(size=n);allocate(s(n));do i=1,n;s(i)=base+79*i;end do;call random_seed(put=s)
   end subroutine seed_rng

   real(dp) function randn() result(z)
      real(dp) :: u1,u2
      call random_number(u1);call random_number(u2);u1=max(u1,1.0e-12_dp)
      z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
   end function randn

   subroutine test_gaussian_copula_gamma
      integer,parameter :: n=64
      real(dp) :: y(n),x(n,2),xs(n,1),time(n),zz,u,mu,rho,sig
      integer :: group(n),i
      type(correlation_spec) :: cor
      type(gaussian_copula_result_t) :: fit
      type(gamlss_control_t) :: ctl
      call seed_rng(601)
      rho=0.58_dp;sig=0.34_dp;zz=randn()
      do i=1,n
         x(i,1)=1.0_dp;x(i,2)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp);xs(i,1)=1.0_dp
         time(i)=real(i,dp);group(i)=1
         if(i>1)zz=rho*zz+sqrt(1.0_dp-rho*rho)*randn()
         u=min(1.0_dp-1.0e-9_dp,max(1.0e-9_dp,normal_cdf(zz)))
         mu=exp(0.35_dp+0.55_dp*x(i,2));y(i)=qGA(u,mu,sig)
      end do
      cor%kind=COR_AR1;cor%fixed=.false.;allocate(cor%par(1));cor%par=0.25_dp
      ctl=gamlss_control_t();ctl%n_cyc=15;ctl%inner_cyc=30;ctl%c_crit=1.0e-4_dp
      call fit_gamlss_gaussian_copula(y,x,GAMLSS_GA,fit,correlation=cor,x_sigma=xs,time=time,group=group, &
         control=ctl,max_iter=100,tolerance=5.0e-6_dp)
      call assert_true(fit%status==0,'Gaussian-copula Gamma status')
      call assert_true(size(fit%correlation_parameters)==1,'Gaussian-copula correlation storage')
      call assert_true(abs(fit%correlation_parameters(1)-rho)<0.25_dp,'Gaussian-copula AR1 recovery')
      call assert_true(abs(fit%model%mu%coefficients(2)-0.55_dp)<0.25_dp,'Gaussian-copula slope recovery')
      call assert_true(fit%copula_log_likelihood>0.0_dp,'Gaussian-copula dependence contribution')
   end subroutine test_gaussian_copula_gamma

   subroutine test_cross_parameter_random_covariance
      integer,parameter :: ng=8,m=13,n=ng*m
      real(dp) :: y(n),xm(n,2),xs(n,1),zr(n,1,2),umu(ng),usig(ng),xx,lsig
      integer :: group(n),g,i,k
      logical :: act(4)
      type(joint_random_effects_result_t) :: fit
      type(gamlss_control_t) :: ctl
      real(dp) :: ce
      call seed_rng(602)
      do g=1,ng
         umu(g)=0.42_dp*sin(0.72_dp*real(g,dp))
         usig(g)=0.45_dp*umu(g)+0.08_dp*cos(0.53_dp*real(g,dp))
      end do
      k=0
      do g=1,ng;do i=1,m
         k=k+1;group(k)=g;xx=-1.0_dp+2.0_dp*real(i-1,dp)/real(m-1,dp)
         xm(k,:)=[1.0_dp,xx];xs(k,1)=1.0_dp;zr(k,1,1)=1.0_dp;zr(k,1,2)=1.0_dp
         lsig=-1.05_dp+usig(g);y(k)=1.20_dp+0.65_dp*xx+umu(g)+exp(lsig)*randn()
      end do;end do
      act=[.true.,.true.,.false.,.false.]
      ctl=gamlss_control_t();ctl%n_cyc=15;ctl%inner_cyc=30;ctl%c_crit=2.0e-4_dp
      call fit_gamlss_joint_random_effects(y,xm,zr,group,GAMLSS_NO,fit,active_parameters=act,x_sigma=xs, &
         control=ctl,max_outer=5,max_inner=70,tol_cov=3.0e-4_dp,tolerance=5.0e-6_dp)
      call assert_true(fit%status==0,'joint random-effects status')
      call assert_true(all(shape(fit%joint_covariance)==[2,2]),'joint random covariance shape')
      call assert_true(fit%joint_covariance(1,2)>0.0_dp,'cross-parameter covariance positive')
      ce=correlation(fit%effects(:,1,1),fit%effects(:,1,2))
      call assert_true(ce>0.20_dp,'fitted mu/sigma random effects correlated')
      call assert_true(abs(fit%model%mu%coefficients(2)-0.65_dp)<0.22_dp,'joint random fixed slope')
   end subroutine test_cross_parameter_random_covariance

   real(dp) function correlation(a,b) result(c)
      real(dp),intent(in) :: a(:),b(:)
      real(dp) :: am,bm,da,db
      am=sum(a)/real(size(a),dp);bm=sum(b)/real(size(b),dp)
      da=sum((a-am)**2);db=sum((b-bm)**2)
      if(da<=0.0_dp.or.db<=0.0_dp)then;c=0.0_dp
      else;c=sum((a-am)*(b-bm))/sqrt(da*db);end if
   end function correlation
end program test_v06
