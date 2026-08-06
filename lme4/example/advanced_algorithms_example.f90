program advanced_algorithms_example
   use lme4
   implicit none

   call demonstrate_pls_lmm()
   call demonstrate_multidimensional_aghq()
   call demonstrate_nlmm()

contains

   subroutine demonstrate_pls_lmm()
      integer, parameter :: ng=6,m=6,n=ng*m
      real(dp) :: y(n),x(n,2),z(n,1),b(ng),t
      type(random_term_t) :: terms(1)
      type(lmm_result_t) :: fit
      integer :: i,g,j

      b=[-0.6_dp,-0.3_dp,-0.1_dp,0.1_dp,0.3_dp,0.6_dp]
      allocate(terms(1)%group(n))
      do g=1,ng
         do j=1,m
            i=(g-1)*m+j
            t=-1.0_dp+2.0_dp*real(j-1,dp)/real(m-1,dp)
            x(i,:)=[1.0_dp,t]
            z(i,1)=1.0_dp
            y(i)=1.2_dp+0.7_dp*t+b(g)+0.05_dp*sin(real(i,dp))
            terms(1)%group(i)=g
         end do
      end do
      terms(1)%z=z
      terms(1)%n_levels=ng
      terms(1)%name='group'
      call fit_lmm_pls(y,x,terms,fit,reml=.false.)
      write(*,'(/,a)') 'Woodbury/penalized least-squares LMM'
      write(*,'(a,2f12.6)') 'fixed effects: ',fit%beta
      write(*,'(a,f12.6)') 'residual sigma: ',fit%sigma
   end subroutine demonstrate_pls_lmm

   subroutine demonstrate_multidimensional_aghq()
      integer, parameter :: ng=6,m=10,n=ng*m
      real(dp) :: x(n,2),z(n,2),beta(2)
      real(dp), allocatable :: y(:),u(:)
      type(random_term_t) :: term,terms(1)
      type(covariance_block_t) :: varcorr(1)
      type(glmm_result_t) :: fit
      type(glmm_control_t) :: control
      integer :: i,g,j

      allocate(term%group(n))
      do g=1,ng
         do j=1,m
            i=(g-1)*m+j
            x(i,:)=[1.0_dp,-1.0_dp+2.0_dp*real(j-1,dp)/real(m-1,dp)]
            z(i,:)=x(i,:)
            term%group(i)=g
         end do
      end do
      term%z=z
      term%n_levels=ng
      term%name='group'
      terms(1)=term
      allocate(varcorr(1)%covariance(2,2),varcorr(1)%sdcor(2,2))
      varcorr(1)%covariance=reshape([0.36_dp,0.08_dp,0.08_dp,0.16_dp],[2,2])
      call cov2sdcor(varcorr(1)%covariance,varcorr(1)%sdcor)
      varcorr(1)%n_levels=ng
      varcorr(1)%name='group'
      beta=[-0.3_dp,0.8_dp]
      call simulate_glmm(x,terms,beta,varcorr,family_binomial,y,u,seed=987)
      control%maxfun=500
      control%tolerance=5.0e-5_dp
      call fit_glmm_aghq_multidimensional(y,x,term,family_binomial,fit,3,control=control)
      write(*,'(/,a)') 'Three-node tensor AGHQ with random intercept and slope'
      write(*,'(a,2f12.6)') 'fixed effects: ',fit%beta
      write(*,'(a,4f12.6)') 'random covariance: ',fit%varcorr(1)%covariance
   end subroutine demonstrate_multidimensional_aghq

   subroutine demonstrate_nlmm()
      integer, parameter :: ng=6,m=7,n=ng*m
      real(dp) :: covariates(n,1),y(n),b(ng),t
      integer :: group(n),i,g,j
      type(nlmm_result_t) :: fit
      type(nlmm_control_t) :: control

      b=[-0.35_dp,-0.20_dp,-0.05_dp,0.05_dp,0.20_dp,0.35_dp]
      do g=1,ng
         do j=1,m
            i=(g-1)*m+j
            t=real(j-1,dp)/real(m-1,dp)
            covariates(i,1)=t
            group(i)=g
            y(i)=1.4_dp*exp(0.55_dp*t)+b(g)+0.025_dp*sin(1.3_dp*real(i,dp))
         end do
      end do
      control%maxfun=700
      control%tolerance=1.0e-5_dp
      call fit_nlmm(y,covariates,group,ng,1,growth_curve,[1.2_dp,0.4_dp],fit, &
         control=control,sigma_start=0.1_dp)
      write(*,'(/,a)') 'Gaussian nonlinear mixed model'
      write(*,'(a,2f12.6)') 'nonlinear fixed effects: ',fit%beta
      write(*,'(a,f12.6)') 'residual sigma: ',fit%sigma
   end subroutine demonstrate_nlmm

   function growth_curve(covariates,beta,random_effect) result(mean)
      real(dp), intent(in) :: covariates(:),beta(:),random_effect(:)
      real(dp) :: mean
      mean=beta(1)*exp(beta(2)*covariates(1))+random_effect(1)
   end function growth_curve

end program advanced_algorithms_example
