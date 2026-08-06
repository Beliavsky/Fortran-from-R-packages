program glmm_extensions_example
   use lme4
   implicit none
   integer, parameter :: ng=8, m=12, n=ng*m
   real(dp) :: x(n,2), z(n,1), beta(2)
   real(dp), allocatable :: y(:), u(:)
   type(random_term_t) :: terms(1)
   type(covariance_block_t) :: vc(1)
   type(glmm_result_t) :: nb_fit, aghq_fit
   type(glmm_control_t) :: control
   integer :: g, j, i

   allocate(terms(1)%group(n))
   do g = 1, ng
      do j = 1, m
         i = (g-1)*m+j
         x(i,:) = [1.0_dp,-1.0_dp+2.0_dp*real(j-1,dp)/real(m-1,dp)]
         z(i,1) = 1.0_dp
         terms(1)%group(i) = g
      end do
   end do
   terms(1)%z = z
   terms(1)%n_levels = ng
   terms(1)%name = 'cluster'

   allocate(vc(1)%covariance(1,1),vc(1)%sdcor(1,1))
   vc(1)%covariance(1,1) = 0.25_dp
   vc(1)%sdcor(1,1) = 0.50_dp
   vc(1)%name = 'cluster'
   vc(1)%n_levels = ng

   beta = [0.25_dp,0.40_dp]
   call simulate_glmm(x,terms,beta,vc,family_negative_binomial,y,u, &
      seed=2718,dispersion=3.0_dp)
   control%maxfun = 400
   control%max_profile = 20
   control%lower_log_dispersion = log(0.25_dp)
   control%upper_log_dispersion = log(30.0_dp)
   call fit_glmer_nb(y,x,terms,nb_fit,control=control)
   if (.not. nb_fit%converged) error stop nb_fit%message
   write(*,'(a)') 'Profiled negative-binomial GLMM'
   write(*,'(a,2f12.6)') 'fixed effects: ',nb_fit%beta
   write(*,'(a,f12.6)') 'random-effect SD: ',nb_fit%theta(1)
   write(*,'(a,f12.6)') 'profiled size: ',nb_fit%dispersion

   deallocate(y,u)
   beta = [-0.30_dp,0.85_dp]
   call simulate_glmm(x,terms,beta,vc,family_binomial,y,u,seed=3141)
   call fit_glmm_aghq(y,x,terms,family_binomial,aghq_fit,order=7,control=control)
   if (.not. aghq_fit%converged) error stop aghq_fit%message
   write(*,'(/,a)') 'Seven-node adaptive-quadrature binomial GLMM'
   write(*,'(a,2f12.6)') 'fixed effects: ',aghq_fit%beta
   write(*,'(a,f12.6)') 'random-effect SD: ',aghq_fit%theta(1)
   write(*,'(a,f12.6)') 'marginal deviance: ',aghq_fit%deviance
end program glmm_extensions_example
