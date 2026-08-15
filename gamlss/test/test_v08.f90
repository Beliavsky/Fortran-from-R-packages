program test_v08
   use gamlss
   use gamlss_kinds, only : dp
   use gamlss_fit, only : GAMLSS_NO
   use gamlss_mvn_v07, only : mvn_logpdf
   implicit none
   integer :: failures
   failures=0
   call test_ais_gaussian_oracle(failures)
   call test_ais_dimension_six(failures)
   call test_ais_refinement(failures)
   if(failures/=0)then
      write(*,'(a,i0)') 'test_v08: FAIL ',failures
      error stop 1
   end if
   write(*,'(a)') 'test_v08: PASS'
contains
   subroutine test_ais_gaussian_oracle(failures)
      integer,intent(inout) :: failures
      integer,parameter :: ng=4,m=5,n=ng*m
      real(dp) :: y(n),xmu(n,2),xs(n,1),z(n,1,2),cov0(1,1)
      integer :: grp(n),g,i,j,status
      real(dp) :: x,mu,sig,b,eps,ll_exact,lp
      real(dp),allocatable :: yy(:),mm(:),cv(:,:)
      logical :: active(4)
      type(joint_random_ais_result_t) :: fit
      j=0
      do g=1,ng
         b=0.18_dp*real(g-2,dp)
         do i=1,m
            j=j+1;x=-0.8_dp+1.6_dp*real(i-1,dp)/real(m-1,dp)
            xmu(j,:)=[1.0_dp,x];xs(j,1)=1.0_dp;grp(j)=g
            z(j,1,1)=1.0_dp;z(j,1,2)=0.0_dp
            mu=0.7_dp+0.55_dp*x+b;sig=0.42_dp
            eps=0.12_dp*sin(real(3*j+g,dp))+0.04_dp*cos(real(2*j,dp))
            y(j)=mu+sig*eps
         end do
      end do
      cov0(1,1)=0.18_dp**2;active=.false.;active(1)=.true.
      call fit_gamlss_joint_random_effects_ais(y,xmu,z,grp,GAMLSS_NO,fit,active_parameters=active, &
         x_sigma=xs,initial_covariance=cov0,qmc_points=2048,proposal_scale=1.15_dp)
      if(fit%status/=0.or..not.(fit%marginal_log_likelihood>-huge(1.0_dp)/100.0_dp))then
         failures=failures+1;return
      end if
      sig=fit%model%sigma%fitted(1);ll_exact=0.0_dp
      allocate(yy(m),mm(m),cv(m,m));cv=fit%joint_covariance(1,1)
      do i=1,m;cv(i,i)=cv(i,i)+sig*sig;end do
      do g=1,ng
         yy=y((g-1)*m+1:g*m)
         mm=matmul(xmu((g-1)*m+1:g*m,:),fit%model%mu%coefficients)
         lp=mvn_logpdf(yy,mm,cv,status)
         if(status/=0)then;failures=failures+1;return;end if
         ll_exact=ll_exact+lp
      end do
      if(abs(fit%marginal_log_likelihood-ll_exact)>2.0e-2_dp)failures=failures+1
      if(fit%minimum_ess<100.0_dp)failures=failures+1
   end subroutine test_ais_gaussian_oracle

   subroutine test_ais_dimension_six(failures)
      integer,intent(inout) :: failures
      integer,parameter :: ng=3,m=4,n=ng*m,q=3
      real(dp) :: y(n),xmu(n,2),xs(n,1),z(n,q,2),cov0(2*q,2*q)
      integer :: grp(n),g,i,j,k
      real(dp) :: x,mu,sig
      logical :: active(4)
      type(joint_random_ais_result_t) :: fit
      j=0;cov0=0.0_dp
      do k=1,2*q;cov0(k,k)=0.025_dp+0.003_dp*real(k,dp);end do
      do g=1,ng
         do i=1,m
            j=j+1;x=-0.7_dp+1.4_dp*real(i-1,dp)/real(m-1,dp)
            grp(j)=g;xmu(j,:)=[1.0_dp,x];xs(j,1)=1.0_dp
            z(j,:,1)=[1.0_dp,x,x*x];z(j,:,2)=[1.0_dp,x,x*x]
            mu=0.5_dp+0.4_dp*x+0.05_dp*real(g-2,dp)*(1.0_dp+x)
            sig=exp(-0.75_dp+0.03_dp*real(g-2,dp))
            y(j)=mu+sig*(0.15_dp*sin(real(j+2*g,dp))+0.05_dp*cos(real(2*j,dp)))
         end do
      end do
      active=.false.;active(1:2)=.true.
      call fit_gamlss_joint_random_effects_ais(y,xmu,z,grp,GAMLSS_NO,fit,active_parameters=active, &
         x_sigma=xs,initial_covariance=cov0,qmc_points=256,proposal_scale=1.35_dp)
      if(fit%status/=0)then;failures=failures+1;return;end if
      if(fit%latent_dimension/=6)failures=failures+1
      if(.not.(fit%marginal_log_likelihood>-huge(1.0_dp)/100.0_dp))failures=failures+1
      if(any(fit%effective_sample_size<=1.0_dp))failures=failures+1
      if(any(shape(fit%posterior_covariance)/=[ng,6,6]))failures=failures+1
   end subroutine test_ais_dimension_six

   subroutine test_ais_refinement(failures)
      integer,intent(inout) :: failures
      real(dp) :: y(8),x(8,1),xs(8,1),z(8,1,2),c(1,1)
      integer :: g(8),i
      logical :: a(4)
      type(joint_random_ais_result_t) :: fit
      x=1.0_dp;xs=1.0_dp;z=0.0_dp;z(:,1,1)=1.0_dp;c(1,1)=0.04_dp
      do i=1,8
         g(i)=1+(i-1)/4
         y(i)=0.8_dp+0.12_dp*real(g(i)-1,dp)+0.2_dp*sin(real(i,dp))
      end do
      a=.false.;a(1)=.true.
      call fit_gamlss_joint_random_effects_ais(y,x,z,g,GAMLSS_NO,fit,active_parameters=a,x_sigma=xs, &
         initial_covariance=c,qmc_points=64,refine_parameters=.true.,max_iter=1)
      if(fit%status/=0.or..not.fit%parameters_refined)failures=failures+1
      if(fit%optimizer_status/=0.and.fit%optimizer_status/=1)failures=failures+1
   end subroutine test_ais_refinement
end program test_v08
