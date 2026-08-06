program test_lme4
   use lme4
   implicit none
   integer :: failures

   failures = 0
   call test_covariance(failures)
   call test_quadrature(failures)
   call test_lmm_random_intercept(failures)
   call test_lmm_random_slope(failures)
   call test_multiple_terms(failures)
   call test_binomial_glmm(failures)
   call test_poisson_glmm(failures)
   call test_gamma_glmm(failures)
   call test_inverse_gaussian_glmm(failures)
   call test_negative_binomial_glmm(failures)
   call test_aghq_glmm(failures)
   call test_structured_covariance(failures)
   call test_lmm_pls(failures)
   call test_custom_family(failures)
   call test_multidimensional_aghq(failures)
   call test_nlmm(failures)
   call test_inference_and_grouped(failures)

   if (failures /= 0) then
      write(*,'(a,i0)') 'FAILED tests: ', failures
      error stop 1
   end if
   write(*,'(a)') 'All lme4-fortran tests passed.'

contains

   subroutine check(condition,name,failures)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      integer, intent(inout) :: failures
      if (.not. condition) then
         failures = failures+1
         write(*,'(a)') 'FAIL: '//trim(name)
      end if
   end subroutine check

   subroutine test_covariance(failures)
      integer, intent(inout) :: failures
      real(dp) :: s(2,2)
      real(dp), allocatable :: v(:,:), s2(:,:)
      s = reshape([2.0_dp,0.25_dp,0.25_dp,3.0_dp],[2,2])
      call sdcor2cov(s,v)
      call cov2sdcor(v,s2)
      call check(maxval(abs(s-s2)) < 1.0e-12_dp,'sdcor covariance round trip',failures)
      call check(abs(v(1,2)-1.5_dp) < 1.0e-12_dp,'sdcor covariance off diagonal',failures)
   end subroutine test_covariance

   subroutine test_quadrature(failures)
      integer, intent(inout) :: failures
      type(gh_rule_t) :: rule
      integer :: info
      call gh_rule(9,rule,info)
      call check(info == 0,'GH rule status',failures)
      call check(abs(sum(rule%weights)-1.0_dp) < 1.0e-12_dp,'GH weights sum',failures)
      call check(abs(sum(rule%weights*rule%nodes)) < 1.0e-12_dp,'GH first moment',failures)
      call check(abs(sum(rule%weights*rule%nodes**2)-1.0_dp) < 1.0e-10_dp,'GH second moment',failures)
   end subroutine test_quadrature

   subroutine test_lmm_random_intercept(failures)
      integer, intent(inout) :: failures
      integer, parameter :: ng=10, m=6, n=ng*m
      real(dp) :: y(n), x(n,2), zr(n,1), b(ng), e
      type(random_term_t) :: terms(1)
      type(lmm_result_t) :: fit
      integer :: i, g, j

      b = [-1.1_dp,-0.8_dp,-0.5_dp,-0.2_dp,0.0_dp,0.1_dp,0.3_dp,0.5_dp,0.8_dp,1.0_dp]
      allocate(terms(1)%group(n))
      do g=1,ng
         do j=1,m
            i=(g-1)*m+j
            x(i,1)=1.0_dp
            x(i,2)=real(j-1,dp)/real(m-1,dp)-0.5_dp
            zr(i,1)=1.0_dp
            e=0.12_dp*sin(1.7_dp*real(i,dp))
            y(i)=1.5_dp+0.8_dp*x(i,2)+b(g)+e
            terms(1)%group(i)=g
         end do
      end do
      terms(1)%z=zr
      terms(1)%n_levels=ng
      terms(1)%name='subject'
      call fit_lmm(y,x,terms,fit,reml=.false.)
      call check(fit%converged,'random-intercept LMM converged',failures)
      if (fit%converged) then
         call check(abs(fit%beta(1)-1.51_dp) < 0.12_dp,'LMM intercept',failures)
         call check(abs(fit%beta(2)-0.8_dp) < 0.08_dp,'LMM slope',failures)
         call check(fit%varcorr(1)%covariance(1,1) > 0.15_dp,'LMM random variance',failures)
         call check(fit%sigma > 0.03_dp .and. fit%sigma < 0.25_dp,'LMM residual sigma',failures)
      end if
   end subroutine test_lmm_random_intercept

   subroutine test_lmm_random_slope(failures)
      integer, intent(inout) :: failures
      integer, parameter :: ng=8, m=9, n=ng*m
      real(dp) :: y(n), x(n,2), zr(n,2), t, e
      real(dp) :: b0(ng), b1(ng)
      type(random_term_t) :: terms(1)
      type(lmm_result_t) :: fit
      type(lmm_control_t) :: ctrl
      integer :: i, g, j

      b0=[-0.8_dp,-0.5_dp,-0.2_dp,0.0_dp,0.1_dp,0.3_dp,0.5_dp,0.7_dp]
      b1=[-0.35_dp,-0.20_dp,-0.08_dp,0.0_dp,0.05_dp,0.12_dp,0.22_dp,0.30_dp]
      allocate(terms(1)%group(n))
      do g=1,ng
         do j=1,m
            i=(g-1)*m+j
            t=-1.0_dp+2.0_dp*real(j-1,dp)/real(m-1,dp)
            x(i,:)=[1.0_dp,t]
            zr(i,:)=[1.0_dp,t]
            e=0.08_dp*cos(1.31_dp*real(i,dp))
            y(i)=2.0_dp+0.55_dp*t+b0(g)+b1(g)*t+e
            terms(1)%group(i)=g
         end do
      end do
      terms(1)%z=zr
      terms(1)%n_levels=ng
      terms(1)%name='subject'
      ctrl%maxfun=1000
      ctrl%tolerance=1.0e-6_dp
      call fit_lmm(y,x,terms,fit,reml=.false.,control=ctrl)
      call check(fit%converged,'random-slope LMM converged',failures)
      if (fit%converged) then
         call check(abs(fit%beta(1)-2.0125_dp) < 0.15_dp,'random-slope intercept',failures)
         call check(abs(fit%beta(2)-0.5575_dp) < 0.12_dp,'random-slope fixed slope',failures)
         call check(all([(fit%varcorr(1)%covariance(i,i)>0.0_dp,i=1,2)]),'positive random-slope variances',failures)
      end if
   end subroutine test_lmm_random_slope


   subroutine test_multiple_terms(failures)
      integer, intent(inout) :: failures
      integer, parameter :: ns=8, nb=5, m=10, n=ns*m
      real(dp) :: y(n), x(n,2), z1(n,1), z2(n,1), t, e
      real(dp) :: subject_effect(ns), batch_effect(nb)
      type(random_term_t) :: terms(2)
      type(lmm_result_t) :: fit
      type(lmm_control_t) :: ctrl
      integer :: i, s, j, batch

      subject_effect=[-0.7_dp,-0.4_dp,-0.2_dp,0.0_dp,0.1_dp,0.25_dp,0.4_dp,0.6_dp]
      batch_effect=[-0.35_dp,-0.15_dp,0.0_dp,0.15_dp,0.35_dp]
      allocate(terms(1)%group(n),terms(2)%group(n))
      do s=1,ns
         do j=1,m
            i=(s-1)*m+j
            batch=mod(i-1,nb)+1
            t=-1.0_dp+2.0_dp*real(j-1,dp)/real(m-1,dp)
            x(i,:)=[1.0_dp,t]
            z1(i,1)=1.0_dp
            z2(i,1)=1.0_dp
            e=0.06_dp*sin(0.91_dp*real(i,dp))
            y(i)=1.25_dp+0.45_dp*t+subject_effect(s)+batch_effect(batch)+e
            terms(1)%group(i)=s
            terms(2)%group(i)=batch
         end do
      end do
      terms(1)%z=z1
      terms(1)%n_levels=ns
      terms(1)%name='subject'
      terms(2)%z=z2
      terms(2)%n_levels=nb
      terms(2)%name='batch'
      ctrl%maxfun=800
      call fit_lmm(y,x,terms,fit,reml=.false.,control=ctrl)
      call check(fit%converged,'multiple-term LMM converged',failures)
      if (fit%converged) then
         call check(size(fit%varcorr)==2,'multiple-term variance blocks',failures)
         call check(abs(fit%beta(2)-0.45_dp)<0.08_dp,'multiple-term fixed slope',failures)
         call check(all([(fit%varcorr(i)%covariance(1,1)>1.0e-4_dp,i=1,2)]), &
            'multiple-term positive variances',failures)
      end if
   end subroutine test_multiple_terms

   subroutine test_binomial_glmm(failures)
      integer, intent(inout) :: failures
      integer, parameter :: ng=12, m=15, n=ng*m
      real(dp) :: x(n,2), zr(n,1), beta(2)
      real(dp), allocatable :: y(:), u(:)
      type(random_term_t) :: terms(1)
      type(covariance_block_t) :: vc(1)
      type(glmm_result_t) :: fit
      type(glmm_control_t) :: ctrl
      integer :: i,g,j

      allocate(terms(1)%group(n))
      do g=1,ng
         do j=1,m
            i=(g-1)*m+j
            x(i,:)=[1.0_dp,-1.2_dp+2.4_dp*real(j-1,dp)/real(m-1,dp)]
            zr(i,1)=1.0_dp
            terms(1)%group(i)=g
         end do
      end do
      terms(1)%z=zr
      terms(1)%n_levels=ng
      terms(1)%name='cluster'
      allocate(vc(1)%covariance(1,1),vc(1)%sdcor(1,1))
      vc(1)%covariance(1,1)=0.49_dp
      vc(1)%sdcor(1,1)=0.7_dp
      vc(1)%name='cluster'
      vc(1)%n_levels=ng
      beta=[-0.35_dp,1.05_dp]
      call simulate_glmm(x,terms,beta,vc,family_binomial,y,u,seed=441)
      ctrl%maxfun=500
      ctrl%max_pirls=80
      call fit_glmm(y,x,terms,family_binomial,fit,control=ctrl)
      call check(fit%converged,'binomial GLMM converged',failures)
      if (fit%converged) then
         call check(abs(fit%beta(2)-beta(2)) < 0.65_dp,'binomial GLMM slope',failures)
         call check(all(fit%fitted>0.0_dp .and. fit%fitted<1.0_dp),'binomial fitted range',failures)
      end if
   end subroutine test_binomial_glmm

   subroutine test_poisson_glmm(failures)
      integer, intent(inout) :: failures
      integer, parameter :: ng=10, m=12, n=ng*m
      real(dp) :: x(n,2), zr(n,1), beta(2)
      real(dp), allocatable :: y(:), u(:)
      type(random_term_t) :: terms(1)
      type(covariance_block_t) :: vc(1)
      type(glmm_result_t) :: fit
      type(glmm_control_t) :: ctrl
      integer :: i,g,j

      allocate(terms(1)%group(n))
      do g=1,ng
         do j=1,m
            i=(g-1)*m+j
            x(i,:)=[1.0_dp,-1.0_dp+2.0_dp*real(j-1,dp)/real(m-1,dp)]
            zr(i,1)=1.0_dp
            terms(1)%group(i)=g
         end do
      end do
      terms(1)%z=zr
      terms(1)%n_levels=ng
      terms(1)%name='site'
      allocate(vc(1)%covariance(1,1),vc(1)%sdcor(1,1))
      vc(1)%covariance(1,1)=0.25_dp
      vc(1)%sdcor(1,1)=0.5_dp
      vc(1)%name='site'
      vc(1)%n_levels=ng
      beta=[0.45_dp,0.38_dp]
      call simulate_glmm(x,terms,beta,vc,family_poisson,y,u,seed=927)
      ctrl%maxfun=500
      ctrl%max_pirls=80
      call fit_glmm(y,x,terms,family_poisson,fit,control=ctrl)
      call check(fit%converged,'poisson GLMM converged',failures)
      if (fit%converged) then
         call check(abs(fit%beta(2)-beta(2)) < 0.45_dp,'poisson GLMM slope',failures)
         call check(all(fit%fitted>0.0_dp),'poisson positive fitted values',failures)
      end if
   end subroutine test_poisson_glmm

   subroutine setup_extended_glmm(terms,vc,x,z,beta,ng,m)
      type(random_term_t), intent(out) :: terms(:)
      type(covariance_block_t), intent(out) :: vc(:)
      real(dp), allocatable, intent(out) :: x(:,:), z(:,:), beta(:)
      integer, intent(in) :: ng, m
      integer :: n, g, j, i

      n = ng*m
      allocate(x(n,2),z(n,1),beta(2))
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
      terms(1)%name = 'group'
      allocate(vc(1)%covariance(1,1),vc(1)%sdcor(1,1))
      vc(1)%covariance(1,1) = 0.16_dp
      vc(1)%sdcor(1,1) = 0.4_dp
      vc(1)%name = 'group'
      vc(1)%n_levels = ng
      beta = [0.35_dp,0.45_dp]
   end subroutine setup_extended_glmm

   subroutine test_gamma_glmm(failures)
      integer, intent(inout) :: failures
      type(random_term_t) :: terms(1)
      type(covariance_block_t) :: vc(1)
      real(dp), allocatable :: x(:,:), z(:,:), beta(:), y(:), u(:)
      type(glmm_result_t) :: fit
      type(glmm_control_t) :: ctrl
      real(dp), parameter :: dispersion = 0.20_dp

      call setup_extended_glmm(terms,vc,x,z,beta,10,12)
      call simulate_glmm(x,terms,beta,vc,family_gamma,y,u,seed=101, &
         dispersion=dispersion)
      ctrl%maxfun = 350
      ctrl%max_pirls = 100
      ctrl%tolerance = 1.0e-5_dp
      call fit_glmm(y,x,terms,family_gamma,fit,control=ctrl,dispersion=dispersion)
      call check(fit%converged,'gamma GLMM converged',failures)
      if (fit%converged) then
         call check(abs(fit%beta(2)-beta(2)) < 0.25_dp,'gamma GLMM slope',failures)
         call check(all(fit%fitted > 0.0_dp),'gamma fitted values positive',failures)
         call check(abs(fit%dispersion-dispersion) < 1.0e-14_dp, &
            'gamma dispersion retained',failures)
      end if
   end subroutine test_gamma_glmm

   subroutine test_inverse_gaussian_glmm(failures)
      integer, intent(inout) :: failures
      type(random_term_t) :: terms(1)
      type(covariance_block_t) :: vc(1)
      real(dp), allocatable :: x(:,:), z(:,:), beta(:), y(:), u(:)
      type(glmm_result_t) :: fit
      type(glmm_control_t) :: ctrl
      real(dp), parameter :: dispersion = 0.15_dp

      call setup_extended_glmm(terms,vc,x,z,beta,10,12)
      call simulate_glmm(x,terms,beta,vc,family_inverse_gaussian,y,u,seed=102, &
         dispersion=dispersion)
      ctrl%maxfun = 350
      ctrl%max_pirls = 100
      ctrl%tolerance = 1.0e-5_dp
      call fit_glmm(y,x,terms,family_inverse_gaussian,fit,control=ctrl, &
         dispersion=dispersion)
      call check(fit%converged,'inverse-Gaussian GLMM converged',failures)
      if (fit%converged) then
         call check(abs(fit%beta(2)-beta(2)) < 0.25_dp, &
            'inverse-Gaussian GLMM slope',failures)
         call check(all(fit%fitted > 0.0_dp), &
            'inverse-Gaussian fitted values positive',failures)
      end if
   end subroutine test_inverse_gaussian_glmm

   subroutine test_negative_binomial_glmm(failures)
      integer, intent(inout) :: failures
      type(random_term_t) :: terms(1)
      type(covariance_block_t) :: vc(1)
      real(dp), allocatable :: x(:,:), z(:,:), beta(:), y(:), u(:)
      type(glmm_result_t) :: fixed_fit, endpoint_fit, profiled_fit
      type(glmm_control_t) :: ctrl
      real(dp), parameter :: true_size = 4.0_dp

      call setup_extended_glmm(terms,vc,x,z,beta,10,15)
      call simulate_glmm(x,terms,beta,vc,family_negative_binomial,y,u,seed=201, &
         dispersion=true_size)
      ctrl%maxfun = 300
      ctrl%max_pirls = 100
      ctrl%tolerance = 2.0e-5_dp
      ctrl%max_profile = 18
      ctrl%lower_log_dispersion = log(0.5_dp)
      ctrl%upper_log_dispersion = log(20.0_dp)
      call fit_glmm(y,x,terms,family_negative_binomial,fixed_fit,control=ctrl, &
         dispersion=true_size)
      call fit_glmm(y,x,terms,family_negative_binomial,endpoint_fit,control=ctrl, &
         dispersion=0.5_dp)
      call fit_glmer_nb(y,x,terms,profiled_fit,control=ctrl)
      call check(fixed_fit%converged,'fixed-size negative-binomial GLMM converged',failures)
      call check(profiled_fit%converged,'profiled negative-binomial GLMM converged',failures)
      if (fixed_fit%converged) then
         call check(abs(fixed_fit%beta(2)-beta(2)) < 0.40_dp, &
            'negative-binomial GLMM slope',failures)
      end if
      if (profiled_fit%converged .and. endpoint_fit%converged) then
         call check(profiled_fit%dispersion >= 0.5_dp .and. &
            profiled_fit%dispersion <= 20.0_dp,'profiled size within bounds',failures)
         call check(profiled_fit%deviance <= endpoint_fit%deviance+1.0e-6_dp, &
            'profiled size improves endpoint objective',failures)
      end if
   end subroutine test_negative_binomial_glmm

   subroutine test_aghq_glmm(failures)
      integer, intent(inout) :: failures
      type(random_term_t) :: terms(1)
      type(covariance_block_t) :: vc(1)
      real(dp), allocatable :: x(:,:), z(:,:), beta(:), y(:), u(:)
      type(glmm_result_t) :: binomial_fit, poisson_fit
      type(glmm_control_t) :: ctrl

      call setup_extended_glmm(terms,vc,x,z,beta,8,18)
      beta = [-0.25_dp,0.90_dp]
      vc(1)%covariance(1,1) = 0.36_dp
      vc(1)%sdcor(1,1) = 0.60_dp
      call simulate_glmm(x,terms,beta,vc,family_binomial,y,u,seed=301)
      ctrl%maxfun = 600
      ctrl%max_pirls = 100
      ctrl%tolerance = 2.0e-5_dp
      call fit_glmm_aghq(y,x,terms,family_binomial,binomial_fit,7,control=ctrl)
      call check(binomial_fit%converged,'binomial AGHQ converged',failures)
      if (binomial_fit%converged) then
         call check(binomial_fit%quadrature_order == 7,'AGHQ order recorded',failures)
         call check(abs(binomial_fit%beta(2)-beta(2)) < 0.60_dp, &
            'binomial AGHQ slope',failures)
         call check(all(binomial_fit%fitted > 0.0_dp .and. &
            binomial_fit%fitted < 1.0_dp),'binomial AGHQ fitted range',failures)
      end if

      deallocate(y,u)
      beta = [0.20_dp,0.35_dp]
      vc(1)%covariance(1,1) = 0.25_dp
      vc(1)%sdcor(1,1) = 0.50_dp
      call simulate_glmm(x,terms,beta,vc,family_poisson,y,u,seed=302)
      ctrl%maxfun = 500
      call fit_glmm_aghq(y,x,terms,family_poisson,poisson_fit,5,control=ctrl)
      call check(poisson_fit%converged,'poisson AGHQ converged',failures)
      if (poisson_fit%converged) then
         call check(poisson_fit%quadrature_order == 5, &
            'poisson AGHQ order recorded',failures)
         call check(all(poisson_fit%fitted > 0.0_dp), &
            'poisson AGHQ fitted values positive',failures)
      end if
   end subroutine test_aghq_glmm


   subroutine test_structured_covariance(failures)
      integer, intent(inout) :: failures
      type(random_term_t) :: term
      real(dp), allocatable :: covariance(:,:)
      integer :: info

      allocate(term%z(2,3),term%group(2))
      term%z = 0.0_dp
      term%group = 1
      term%n_levels = 1

      term%covariance_structure = covariance_diagonal
      call term_covariance_from_eta(term,log([1.0_dp,2.0_dp,3.0_dp]),covariance,info)
      call check(info == 0,'diagonal covariance status',failures)
      if (info == 0) then
         call check(maxval(abs(diagonal(covariance)-[1.0_dp,4.0_dp,9.0_dp])) < 1.0e-12_dp, &
            'diagonal covariance values',failures)
         call check(abs(covariance(1,2))+abs(covariance(1,3))+abs(covariance(2,3)) < 1.0e-12_dp, &
            'diagonal covariance off diagonal',failures)
      end if

      term%covariance_structure = covariance_compound_symmetry
      call term_covariance_from_eta(term,[log(2.0_dp),0.0_dp],covariance,info)
      call check(info == 0,'compound-symmetry covariance status',failures)
      if (info == 0) then
         call check(maxval(abs(diagonal(covariance)-4.0_dp)) < 1.0e-12_dp, &
            'compound-symmetry variances',failures)
         call check(maxval(abs(covariance-diag_matrix(diagonal(covariance)))) < 1.0e-12_dp, &
            'compound-symmetry zero-correlation center',failures)
      end if

      term%covariance_structure = covariance_ar1
      call term_covariance_from_eta(term,[log(2.0_dp),atanh(0.5_dp)],covariance,info)
      call check(info == 0,'AR1 covariance status',failures)
      if (info == 0) then
         call check(abs(covariance(1,1)-4.0_dp) < 1.0e-12_dp, &
            'AR1 variance',failures)
         call check(abs(covariance(1,2)-2.0_dp) < 1.0e-12_dp, &
            'AR1 lag-one covariance',failures)
         call check(abs(covariance(1,3)-1.0_dp) < 1.0e-12_dp, &
            'AR1 lag-two covariance',failures)
      end if
   end subroutine test_structured_covariance

   subroutine test_lmm_pls(failures)
      integer, intent(inout) :: failures
      integer, parameter :: ng=6,m=6,n=ng*m
      real(dp) :: y(n),x(n,2),z(n,1),b(ng),t,e
      type(random_term_t) :: dense_terms(1),pls_terms(1)
      type(lmm_result_t) :: dense_fit,pls_fit
      type(lmm_control_t) :: ctrl
      integer :: i,g,j

      b=[-0.6_dp,-0.3_dp,-0.1_dp,0.1_dp,0.3_dp,0.6_dp]
      allocate(dense_terms(1)%group(n))
      do g=1,ng
         do j=1,m
            i=(g-1)*m+j
            t=-1.0_dp+2.0_dp*real(j-1,dp)/real(m-1,dp)
            x(i,:)=[1.0_dp,t]
            z(i,1)=1.0_dp
            e=0.05_dp*sin(real(i,dp))
            y(i)=1.2_dp+0.7_dp*t+b(g)+e
            dense_terms(1)%group(i)=g
         end do
      end do
      dense_terms(1)%z=z
      dense_terms(1)%n_levels=ng
      dense_terms(1)%name='group'
      pls_terms=dense_terms
      ctrl%maxfun=500
      call fit_lmm(y,x,dense_terms,dense_fit,reml=.false.,control=ctrl)
      call fit_lmm_pls(y,x,pls_terms,pls_fit,reml=.false.,control=ctrl)
      call check(dense_fit%converged .and. pls_fit%converged, &
         'dense and PLS LMM convergence',failures)
      if (dense_fit%converged .and. pls_fit%converged) then
         call check(maxval(abs(dense_fit%beta-pls_fit%beta)) < 1.0e-6_dp, &
            'PLS fixed effects match dense fit',failures)
         call check(abs(dense_fit%deviance-pls_fit%deviance) < 1.0e-5_dp, &
            'PLS objective matches dense fit',failures)
         call check(abs(dense_fit%sigma-pls_fit%sigma) < 1.0e-6_dp, &
            'PLS residual scale matches dense fit',failures)
      end if
   end subroutine test_lmm_pls

   subroutine test_custom_family(failures)
      integer, intent(inout) :: failures
      integer, parameter :: ng=10,m=20,n=ng*m
      real(dp) :: x(n,2),z(n,1),y(n),r(n),t,probability,b(ng)
      type(random_term_t) :: terms(1)
      type(glmm_result_t) :: fit
      type(glmm_control_t) :: ctrl
      type(family_spec_t) :: family
      integer :: i,g,j

      b=[-0.75_dp,-0.55_dp,-0.35_dp,-0.15_dp,-0.05_dp,0.10_dp,0.25_dp,0.40_dp,0.55_dp,0.75_dp]
      allocate(terms(1)%group(n))
      call set_random_seed(2718)
      call random_number(r)
      do g=1,ng
         do j=1,m
            i=(g-1)*m+j
            t=-1.25_dp+2.5_dp*real(j-1,dp)/real(m-1,dp)
            x(i,:)=[1.0_dp,t]
            z(i,1)=1.0_dp
            probability=normal_cdf_scalar(-0.25_dp+0.8_dp*t+b(g))
            y(i)=merge(1.0_dp,0.0_dp,r(i)<probability)
            terms(1)%group(i)=g
         end do
      end do
      terms(1)%z=z
      terms(1)%n_levels=ng
      terms(1)%name='cluster'
      family=binomial_probit_family()
      ctrl%maxfun=600
      ctrl%max_pirls=100
      call fit_glmm_custom(y,x,terms,family,fit,control=ctrl)
      call check(fit%converged,'custom probit GLMM converged',failures)
      if (fit%converged) then
         call check(fit%beta(2)>0.2_dp .and. fit%beta(2)<1.6_dp, &
            'custom probit slope direction and scale',failures)
         call check(all(fit%fitted>0.0_dp .and. fit%fitted<1.0_dp), &
            'custom probit fitted range',failures)
      end if
   end subroutine test_custom_family

   subroutine test_multidimensional_aghq(failures)
      integer, intent(inout) :: failures
      integer, parameter :: ng=6,m=10,n=ng*m
      real(dp) :: x(n,2),z(n,2),beta(2)
      real(dp), allocatable :: y(:),u(:)
      type(random_term_t) :: term,terms(1)
      type(covariance_block_t) :: varcorr(1)
      type(glmm_result_t) :: fit
      type(glmm_control_t) :: ctrl
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
      ctrl%maxfun=500
      ctrl%tolerance=5.0e-5_dp
      call fit_glmm_aghq_multidimensional(y,x,term,family_binomial,fit,3,control=ctrl)
      call check(fit%converged,'multidimensional AGHQ converged',failures)
      if (fit%converged) then
         call check(fit%quadrature_order==3,'multidimensional AGHQ order recorded',failures)
         call check(size(fit%varcorr(1)%covariance,1)==2, &
            'multidimensional AGHQ covariance dimension',failures)
         call check(all(diagonal(fit%varcorr(1)%covariance)>0.0_dp), &
            'multidimensional AGHQ positive variances',failures)
         call check(all(fit%fitted>0.0_dp .and. fit%fitted<1.0_dp), &
            'multidimensional AGHQ fitted range',failures)
      end if
   end subroutine test_multidimensional_aghq

   subroutine test_nlmm(failures)
      integer, intent(inout) :: failures
      integer, parameter :: ng=6,m=7,n=ng*m
      real(dp) :: covariates(n,1),y(n),b(ng),t
      integer :: group(n),i,g,j
      type(nlmm_result_t) :: fit
      type(nlmm_control_t) :: ctrl

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
      ctrl%maxfun=700
      ctrl%tolerance=1.0e-5_dp
      ctrl%max_mode_iterations=40
      call fit_nlmm(y,covariates,group,ng,1,nlmm_test_mean,[1.2_dp,0.4_dp],fit, &
         control=ctrl,sigma_start=0.1_dp)
      call check(fit%converged,'nonlinear mixed model converged',failures)
      if (fit%converged) then
         call check(abs(fit%beta(1)-1.4_dp)<0.12_dp,'NLMM amplitude',failures)
         call check(abs(fit%beta(2)-0.55_dp)<0.12_dp,'NLMM growth rate',failures)
         call check(fit%sigma>0.0_dp .and. fit%covariance(1,1)>0.0_dp, &
            'NLMM positive variances',failures)
      end if
   end subroutine test_nlmm

   subroutine test_inference_and_grouped(failures)
      integer, intent(inout) :: failures
      integer, parameter :: ng=5,m=6,n=ng*m
      real(dp) :: y(n),x(n,2),z(n,1),t,b(ng),statistic,p_value
      integer :: group(n),i,g,j
      real(dp), allocatable :: lower(:),upper(:),prediction(:)
      type(random_term_t) :: terms(1)
      type(lmm_result_t) :: fit
      type(profile_result_t) :: wald,profile
      type(bootstrap_result_t) :: bootstrap
      type(influence_result_t) :: influence
      type(lm_list_result_t) :: grouped

      b=[-0.4_dp,-0.2_dp,0.0_dp,0.2_dp,0.4_dp]
      allocate(terms(1)%group(n))
      do g=1,ng
         do j=1,m
            i=(g-1)*m+j
            t=-1.0_dp+2.0_dp*real(j-1,dp)/real(m-1,dp)
            x(i,:)=[1.0_dp,t]
            z(i,1)=1.0_dp
            group(i)=g
            y(i)=1.0_dp+0.6_dp*t+b(g)+0.08_dp*sin(real(i,dp))
            terms(1)%group(i)=g
         end do
      end do
      terms(1)%z=z
      terms(1)%n_levels=ng
      terms(1)%name='group'
      call fit_lmm(y,x,terms,fit,reml=.false.)
      call check(fit%converged,'inference source LMM converged',failures)
      if (.not. fit%converged) return

      call wald_confint_lmm(fit,0.95_dp,wald)
      call check(all(wald%lower<fit%beta .and. fit%beta<wald%upper), &
         'Wald intervals contain estimates',failures)
      call profile_confint_lmm_beta(y,x,terms,fit,0.80_dp,profile)
      call check(profile%status==0 .and. &
         all(profile%lower<profile%estimate .and. profile%estimate<profile%upper), &
         'profile-likelihood intervals contain estimates',failures)
      call parametric_bootstrap_lmm(x,terms,fit,3,bootstrap,seed=22)
      call check(bootstrap%successful>=2,'parametric bootstrap refits',failures)
      call bootstrap_percentile_confint(bootstrap%beta,0.80_dp,lower,upper)
      call check(all(lower<=upper),'bootstrap percentile interval ordering',failures)
      call influence_lmm_groups(y,x,terms,fit,1,influence)
      call check(count(influence%converged)>=ng-1,'group-deletion influence refits',failures)
      call likelihood_ratio_test(fit%deviance+4.0_dp,fit%deviance,1,statistic,p_value)
      call check(abs(statistic-4.0_dp)<1.0e-12_dp,'likelihood-ratio statistic',failures)
      call check(p_value>0.04_dp .and. p_value<0.05_dp,'likelihood-ratio probability',failures)
      call fit_lm_list(y,x,group,ng,grouped)
      call check(all(grouped%converged),'lmList grouped regressions',failures)
      call predict_lm_list(grouped,x,group,prediction)
      call check(size(prediction)==n .and. all(abs(prediction)<huge(1.0_dp)), &
         'lmList prediction dimensions',failures)
   end subroutine test_inference_and_grouped

   pure function diagonal(a) result(d)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: d(min(size(a,1),size(a,2)))
      integer :: i
      do i=1,size(d)
         d(i)=a(i,i)
      end do
   end function diagonal

   pure function diag_matrix(d) result(a)
      real(dp), intent(in) :: d(:)
      real(dp) :: a(size(d),size(d))
      integer :: i
      a=0.0_dp
      do i=1,size(d)
         a(i,i)=d(i)
      end do
   end function diag_matrix

   function nlmm_test_mean(covariates,beta,random_effect) result(mean)
      real(dp), intent(in) :: covariates(:),beta(:),random_effect(:)
      real(dp) :: mean
      mean=beta(1)*exp(beta(2)*covariates(1))+random_effect(1)
   end function nlmm_test_mean

end program test_lme4
