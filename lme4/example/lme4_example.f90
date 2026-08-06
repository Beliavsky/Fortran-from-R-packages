program lme4_example
   use lme4
   implicit none
   integer, parameter :: ng=18, m=10, n=ng*m
   real(dp) :: x(n,2), zr(n,2), beta_true(2)
   real(dp), allocatable :: y(:), u(:)
   type(random_term_t) :: terms(1)
   type(covariance_block_t) :: vc(1)
   type(lmm_result_t) :: fit
   integer :: i,g,j
   real(dp) :: t

   allocate(terms(1)%group(n))
   do g=1,ng
      do j=1,m
         i=(g-1)*m+j
         t=real(j-1,dp)/real(m-1,dp)
         x(i,:)=[1.0_dp,t]
         zr(i,:)=[1.0_dp,t]
         terms(1)%group(i)=g
      end do
   end do
   terms(1)%z=zr
   terms(1)%n_levels=ng
   terms(1)%name='subject'

   allocate(vc(1)%covariance(2,2),vc(1)%sdcor(2,2))
   vc(1)%covariance=reshape([0.64_dp,0.12_dp,0.12_dp,0.16_dp],[2,2])
   call cov2sdcor(vc(1)%covariance,vc(1)%sdcor)
   vc(1)%name='subject'
   vc(1)%n_levels=ng
   beta_true=[2.0_dp,0.75_dp]

   call simulate_lmm(x,terms,beta_true,vc,0.25_dp,y,u,seed=808)
   call fit_lmm(y,x,terms,fit,reml=.true.)
   if (.not. fit%converged) error stop fit%message

   write(*,'(a)') 'Dense linear mixed model'
   write(*,'(a,2f12.6)') 'fixed effects: ',fit%beta
   write(*,'(a,f12.6)') 'residual sigma: ',fit%sigma
   write(*,'(a)') 'random-effect sd/correlation matrix:'
   do i=1,size(fit%varcorr(1)%sdcor,1)
      write(*,'(*(f12.6,1x))') fit%varcorr(1)%sdcor(i,:)
   end do
   write(*,'(a,f12.4)') 'REML criterion: ',fit%deviance
end program lme4_example
