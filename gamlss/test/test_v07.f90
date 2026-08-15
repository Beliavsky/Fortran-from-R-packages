program test_v07
   use gamlss
   use gamlss_discrete, only : qNBI
   use gamlss_special, only : normal_cdf
   use nlme_types, only : correlation_spec,COR_AR1
   implicit none
   call test_mvn_rectangle_oracle
   call test_mixed_atom_support
   call test_discrete_gaussian_copula
   call test_joint_random_ghq
   print *, 'test_v07: PASS'
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
      call random_seed(size=n);allocate(s(n));do i=1,n;s(i)=base+37*i;end do;call random_seed(put=s)
   end subroutine seed_rng

   real(dp) function randn() result(z)
      real(dp) :: u1,u2
      call random_number(u1);call random_number(u2);u1=max(u1,1.0e-12_dp)
      z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
   end function randn

   subroutine test_mvn_rectangle_oracle
      real(dp) :: lo(2),hi(2),mu(2),cov(2,2),p,ref,err,rho
      integer :: status
      rho=0.60_dp;lo=-huge(1.0_dp);hi=0.0_dp;mu=0.0_dp
      cov=reshape([1.0_dp,rho,rho,1.0_dp],[2,2])
      p=mvn_rectangle_probability(lo,hi,mu,cov,status,n_qmc=4096,error_estimate=err)
      ref=0.25_dp+asin(rho)/(2.0_dp*acos(-1.0_dp))
      call assert_true(status==0,'MVN rectangle status')
      call assert_true(abs(p-ref)<2.0e-5_dp,'MVN rectangle bivariate oracle')
      call assert_true(err<1.0e-3_dp,'MVN rectangle error estimate')
   end subroutine test_mvn_rectangle_oracle

   subroutine test_mixed_atom_support
      real(dp) :: left1,p0
      call assert_true(family_observation_is_atom(GAMLSS_BEINF,0.0_dp),'BEINF zero atom')
      call assert_true(family_observation_is_atom(GAMLSS_BEINF,1.0_dp),'BEINF one atom')
      call assert_true(.not.family_observation_is_atom(GAMLSS_BEINF,0.4_dp),'BEINF interior continuous')
      p0=family_cdf(GAMLSS_BEINF,0.0_dp,0.45_dp,0.30_dp,0.25_dp,0.40_dp)
      left1=family_cdf_left(GAMLSS_BEINF,1.0_dp,0.45_dp,0.30_dp,0.25_dp,0.40_dp)
      call assert_true(abs(p0-0.25_dp/1.65_dp)<1.0e-12_dp,'BEINF zero interval mass')
      call assert_true(abs(left1-1.25_dp/1.65_dp)<1.0e-12_dp,'BEINF one left limit')
   end subroutine test_mixed_atom_support

   subroutine test_discrete_gaussian_copula
      integer,parameter :: ng=10,m=3,n=ng*m
      real(dp) :: y(n),xm(n,1),xs(n,1),time(n),z,rho,u
      integer :: group(n),g,i,k
      type(correlation_spec) :: cor
      type(gaussian_copula_mixed_result_t) :: fit
      type(gamlss_control_t) :: ctl
      call seed_rng(711);rho=0.50_dp;k=0
      do g=1,ng
         z=randn()
         do i=1,m
            k=k+1;if(i>1)z=rho*z+sqrt(1.0_dp-rho*rho)*randn()
            u=min(1.0_dp-1.0e-10_dp,max(1.0e-10_dp,normal_cdf(z)))
            y(k)=real(qNBI(u,2.7_dp,0.55_dp),dp)
            xm(k,1)=1.0_dp;xs(k,1)=1.0_dp;time(k)=real(i,dp);group(k)=g
         end do
      end do
      cor%kind=COR_AR1;cor%fixed=.false.;allocate(cor%par(1));cor%par=0.20_dp
      ctl=gamlss_control_t();ctl%n_cyc=12;ctl%inner_cyc=25
      call fit_gamlss_gaussian_copula_mixed(y,xm,GAMLSS_NBI,fit,correlation=cor,x_sigma=xs,time=time, &
         group=group,control=ctl,n_qmc=512,max_iter=60,tolerance=2.0e-5_dp)
      call assert_true(fit%status==0,'discrete Gaussian-copula status')
      call assert_true(size(fit%correlation_parameters)==1,'discrete copula correlation storage')
      call assert_true(abs(fit%correlation_parameters(1)-rho)<0.25_dp,'discrete copula AR1 recovery')
      call assert_true(all(fit%group_log_likelihood<0.0_dp),'discrete copula group likelihoods')
   end subroutine test_discrete_gaussian_copula

   subroutine test_joint_random_ghq
      integer,parameter :: ng=5,m=8,n=ng*m
      real(dp) :: y(n),xm(n,2),xs(n,1),zr(n,1,2),init(2,2),bmu,bsig,xx
      integer :: group(n),g,i,k
      logical :: act(4)
      type(joint_random_ghq_result_t) :: fit
      type(gamlss_control_t) :: ctl
      call seed_rng(702);k=0
      do g=1,ng
         bmu=0.30_dp*sin(0.73_dp*real(g,dp))
         bsig=0.40_dp*bmu+0.04_dp*cos(0.51_dp*real(g,dp))
         do i=1,m
            k=k+1;xx=-1.0_dp+2.0_dp*real(i-1,dp)/real(m-1,dp);group(k)=g
            xm(k,:)=[1.0_dp,xx];xs(k,1)=1.0_dp;zr(k,1,1)=1.0_dp;zr(k,1,2)=1.0_dp
            y(k)=1.10_dp+0.50_dp*xx+bmu+exp(-1.0_dp+bsig)*randn()
         end do
      end do
      act=[.true.,.true.,.false.,.false.]
      init=reshape([0.08_dp,0.02_dp,0.02_dp,0.03_dp],[2,2])
      ctl=gamlss_control_t();ctl%n_cyc=12;ctl%inner_cyc=25
      call fit_gamlss_joint_random_effects_ghq(y,xm,zr,group,GAMLSS_NO,fit,active_parameters=act, &
         x_sigma=xs,initial_covariance=init,control=ctl,quadrature_order=5,max_iter=60,tolerance=1.0e-5_dp)
      call assert_true(fit%status==0,'joint GHQ status')
      call assert_true(fit%latent_dimension==2.and.fit%quadrature_order==5,'joint GHQ dimensions')
      call assert_true(all(shape(fit%joint_covariance)==[2,2]),'joint GHQ covariance shape')
      call assert_true(fit%joint_covariance(1,1)>0.0_dp.and.fit%joint_covariance(2,2)>0.0_dp, &
         'joint GHQ positive variances')
      call assert_true(determinant2(fit%joint_covariance)>0.0_dp,'joint GHQ covariance SPD')
      call assert_true(abs(fit%model%mu%coefficients(2)-0.50_dp)<0.25_dp,'joint GHQ fixed slope')
      call assert_true(all(shape(fit%posterior_effects)==[ng,1,4]),'joint GHQ posterior effects')
   end subroutine test_joint_random_ghq

   real(dp) function determinant2(a) result(d)
      real(dp),intent(in) :: a(:,:)
      d=a(1,1)*a(2,2)-a(1,2)*a(2,1)
   end function determinant2
end program test_v07
