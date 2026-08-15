program v08_extended
   use gamlss
   use gamlss_kinds, only : dp
   use gamlss_fit, only : GAMLSS_NO
   implicit none
   integer,parameter :: ng=3,m=4,n=ng*m,q=3
   real(dp) :: y(n),xmu(n,2),xs(n,1),z(n,q,2),cov0(2*q,2*q)
   integer :: grp(n),g,i,j,k
   real(dp) :: x,mu,sig
   logical :: active(4)
   type(joint_random_ais_result_t) :: fit
   cov0=0.0_dp
   do k=1,2*q;cov0(k,k)=0.03_dp;end do
   j=0
   do g=1,ng
      do i=1,m
         j=j+1;x=-0.7_dp+1.4_dp*real(i-1,dp)/real(m-1,dp)
         grp(j)=g;xmu(j,:)=[1.0_dp,x];xs(j,1)=1.0_dp
         z(j,:,1)=[1.0_dp,x,x*x];z(j,:,2)=[1.0_dp,x,x*x]
         mu=0.6_dp+0.45_dp*x+0.06_dp*real(g-2,dp)*(1.0_dp+x)
         sig=exp(-0.7_dp+0.04_dp*real(g-2,dp))
         y(j)=mu+sig*(0.14_dp*sin(real(j+g,dp))+0.04_dp*cos(real(2*j,dp)))
      end do
   end do
   active=.false.;active(1:2)=.true.
   call fit_gamlss_joint_random_effects_ais(y,xmu,z,grp,GAMLSS_NO,fit,active_parameters=active, &
      x_sigma=xs,initial_covariance=cov0,qmc_points=512,proposal_scale=1.35_dp)
   write(*,'(a,i0)') 'AIS latent dimension: ',fit%latent_dimension
   write(*,'(a,f12.4)') 'AIS marginal log likelihood: ',fit%marginal_log_likelihood
   write(*,'(a,f10.1)') 'Minimum group ESS: ',fit%minimum_ess
end program v08_extended
