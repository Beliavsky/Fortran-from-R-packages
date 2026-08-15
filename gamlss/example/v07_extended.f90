program v07_extended
   use gamlss
   use gamlss_discrete, only : qNBI
   use gamlss_special, only : normal_cdf
   use nlme_types, only : correlation_spec,COR_AR1
   implicit none
   call discrete_copula_example
   call ghq_random_effects_example
contains
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

   subroutine discrete_copula_example
      integer,parameter :: ng=10,m=3,n=ng*m
      real(dp) :: y(n),xm(n,1),xs(n,1),time(n),z,u,rho
      integer :: group(n),g,i,k
      type(correlation_spec) :: cor
      type(gaussian_copula_mixed_result_t) :: fit
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
      call fit_gamlss_gaussian_copula_mixed(y,xm,GAMLSS_NBI,fit,correlation=cor,x_sigma=xs, &
         time=time,group=group,n_qmc=512,max_iter=60,tolerance=2.0e-5_dp)
      write(*,'(a,f9.4)')'Discrete-copula fitted AR(1) rho: ',fit%correlation_parameters(1)
   end subroutine discrete_copula_example

   subroutine ghq_random_effects_example
      integer,parameter :: ng=5,m=8,n=ng*m
      real(dp) :: y(n),xm(n,2),xs(n,1),zr(n,1,2),init(2,2),bmu,bsig,xx
      integer :: group(n),g,i,k
      logical :: active(4)
      type(joint_random_ghq_result_t) :: fit
      call seed_rng(712);k=0
      do g=1,ng
         bmu=0.30_dp*sin(0.73_dp*real(g,dp));bsig=0.40_dp*bmu+0.04_dp*cos(0.51_dp*real(g,dp))
         do i=1,m
            k=k+1;xx=-1.0_dp+2.0_dp*real(i-1,dp)/real(m-1,dp);group(k)=g
            xm(k,:)=[1.0_dp,xx];xs(k,1)=1.0_dp;zr(k,1,1)=1.0_dp;zr(k,1,2)=1.0_dp
            y(k)=1.10_dp+0.50_dp*xx+bmu+exp(-1.0_dp+bsig)*randn()
         end do
      end do
      active=[.true.,.true.,.false.,.false.]
      init=reshape([0.08_dp,0.02_dp,0.02_dp,0.03_dp],[2,2])
      call fit_gamlss_joint_random_effects_ghq(y,xm,zr,group,GAMLSS_NO,fit,active_parameters=active, &
         x_sigma=xs,initial_covariance=init,quadrature_order=5,max_iter=60,tolerance=1.0e-5_dp)
      write(*,'(a,f10.4)')'GHQ marginal log likelihood: ',fit%marginal_log_likelihood
      write(*,'(a,3f10.5)')'GHQ covariance (v11,v12,v22): ',fit%joint_covariance(1,1), &
         fit%joint_covariance(1,2),fit%joint_covariance(2,2)
   end subroutine ghq_random_effects_example
end program v07_extended
