! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch_copula
   use rmgarch_kinds, only : dp
   use rmgarch_math, only : covariance_matrix, normalize_covariance, make_positive_definite, &
      normal_cdf, normal_quantile
   use rmgarch_distributions, only : multivariate_normal_logpdf, multivariate_student_logpdf, &
      standardized_student_pdf, standardized_student_cdf, standardized_student_quantile, &
      standardized_laplace_cdf, clamp_probabilities, &
      random_multivariate_normal, random_multivariate_student
   use rmgarch_optimizer, only : optimizer_result, nelder_mead
   use rmgarch_types, only : dcc_spec, copula_fit_result, &
      dist_gaussian, dist_student, dist_laplace
   use rmgarch_dcc, only : dcc_filter, make_dcc_spec, simulate_dcc
   implicit none
   private

   type :: copula_context
      real(dp), allocatable :: z(:,:)
      integer :: p = 1
      integer :: q = 1
      integer :: g = 0
      integer :: distribution = dist_gaussian
      real(dp) :: shape = 8.0_dp
      logical :: estimate_shape = .false.
      logical :: time_varying = .true.
   end type copula_context

   public :: gaussian_copula_log_density, student_copula_log_density
   public :: copula_log_density, static_copula_log_likelihood
   public :: dynamic_gaussian_copula_log_likelihood
   public :: dynamic_student_copula_log_likelihood
   public :: fit_copula, simulate_static_copula, simulate_dynamic_copula
   public :: parametric_uniform_transform, copula_score_transform
   public :: score_to_uniform_transform

contains

   function gaussian_copula_log_density(z, r, valid) result(value)
      real(dp), intent(in) :: z(:), r(:,:)
      logical, intent(out), optional :: valid
      real(dp) :: value, zero(size(z)), identity(size(z),size(z))
      logical :: ok1, ok2
      integer :: i

      zero = 0.0_dp
      identity = 0.0_dp
      do i = 1, size(z)
         identity(i,i) = 1.0_dp
      end do
      value = multivariate_normal_logpdf(z,zero,r,ok1)- &
         multivariate_normal_logpdf(z,zero,identity,ok2)
      if (.not. (ok1 .and. ok2)) value = -huge(1.0_dp)
      if (present(valid)) valid = ok1 .and. ok2
   end function gaussian_copula_log_density

   function student_copula_log_density(z, r, nu, valid) result(value)
      real(dp), intent(in) :: z(:), r(:,:), nu
      logical, intent(out), optional :: valid
      real(dp) :: value, zero(size(z)), margins
      logical :: ok
      integer :: j

      zero = 0.0_dp
      value = multivariate_student_logpdf(z,zero,r,nu,ok)
      if (ok) then
         margins = 0.0_dp
         do j = 1, size(z)
            margins = margins+log(max(standardized_student_pdf(z(j),nu),tiny(1.0_dp)))
         end do
         value = value-margins
      else
         value = -huge(1.0_dp)
      end if
      if (present(valid)) valid = ok
   end function student_copula_log_density

   function copula_log_density(distribution, z, r, shape, valid) result(value)
      integer, intent(in) :: distribution
      real(dp), intent(in) :: z(:), r(:,:)
      real(dp), intent(in), optional :: shape
      logical, intent(out), optional :: valid
      real(dp) :: value, nu
      logical :: ok

      nu = 8.0_dp
      if (present(shape)) nu = shape
      select case (distribution)
      case (dist_gaussian)
         value = gaussian_copula_log_density(z,r,ok)
      case (dist_student)
         value = student_copula_log_density(z,r,nu,ok)
      case default
         value = -huge(1.0_dp)
         ok = .false.
      end select
      if (present(valid)) valid = ok
   end function copula_log_density

   function static_copula_log_likelihood(z, r, distribution, shape, valid) result(value)
      real(dp), intent(in) :: z(:,:), r(:,:)
      integer, intent(in) :: distribution
      real(dp), intent(in), optional :: shape
      logical, intent(out), optional :: valid
      real(dp) :: value, nu
      integer :: t
      logical :: ok, step_ok

      nu = 8.0_dp
      if (present(shape)) nu = shape
      ok = size(r,1) == size(z,2) .and. size(r,2) == size(z,2)
      value = 0.0_dp
      if (ok) then
         do t = 1, size(z,1)
            value = value+copula_log_density(distribution,z(t,:),r,nu,step_ok)
            ok = ok .and. step_ok
         end do
      end if
      if (.not. ok) value = -huge(1.0_dp)
      if (present(valid)) valid = ok
   end function static_copula_log_likelihood

   function dynamic_gaussian_copula_log_likelihood(z, spec, valid) result(value)
      real(dp), intent(in) :: z(:,:)
      type(dcc_spec), intent(in) :: spec
      logical, intent(out), optional :: valid
      real(dp) :: value
      real(dp), allocatable :: q(:,:,:), r(:,:,:), ll_dummy(:)
      integer :: t
      logical :: ok, step_ok

      allocate(q(size(z,2),size(z,2),size(z,1)),r(size(z,2),size(z,2),size(z,1)), &
         ll_dummy(size(z,1)))
      call dcc_filter(z,spec,q,r,ll_dummy,valid=ok)
      value = 0.0_dp
      if (ok) then
         do t = 1, size(z,1)
            value = value+gaussian_copula_log_density(z(t,:),r(:,:,t),step_ok)
            ok = ok .and. step_ok
         end do
      end if
      if (.not. ok) value = -huge(1.0_dp)
      if (present(valid)) valid = ok
   end function dynamic_gaussian_copula_log_likelihood

   function dynamic_student_copula_log_likelihood(z, spec, nu, valid) result(value)
      real(dp), intent(in) :: z(:,:), nu
      type(dcc_spec), intent(in) :: spec
      logical, intent(out), optional :: valid
      real(dp) :: value
      real(dp), allocatable :: q(:,:,:), r(:,:,:), ll_dummy(:)
      integer :: t
      logical :: ok, step_ok

      allocate(q(size(z,2),size(z,2),size(z,1)),r(size(z,2),size(z,2),size(z,1)), &
         ll_dummy(size(z,1)))
      call dcc_filter(z,spec,q,r,ll_dummy,valid=ok)
      value = 0.0_dp
      if (ok) then
         do t = 1, size(z,1)
            value = value+student_copula_log_density(z(t,:),r(:,:,t),nu,step_ok)
            ok = ok .and. step_ok
         end do
      end if
      if (.not. ok) value = -huge(1.0_dp)
      if (present(valid)) valid = ok
   end function dynamic_student_copula_log_likelihood

   function fit_copula(z, distribution, time_varying, p, q, g, shape, estimate_shape, &
      max_iterations) result(fit)
      real(dp), intent(in) :: z(:,:)
      integer, intent(in), optional :: distribution, p, q, g, max_iterations
      logical, intent(in), optional :: time_varying, estimate_shape
      real(dp), intent(in), optional :: shape
      type(copula_fit_result) :: fit
      type(copula_context) :: context
      type(optimizer_result) :: opt
      real(dp), allocatable :: x0(:), alpha(:), beta(:), gamma(:)
      real(dp) :: fitted_shape, remaining
      integer :: ncoef, maxit, k
      logical :: ok

      if (present(distribution)) context%distribution = distribution
      if (present(time_varying)) context%time_varying = time_varying
      if (present(p)) context%p = max(0,p)
      if (present(q)) context%q = max(0,q)
      if (present(g)) context%g = max(0,g)
      if (present(shape)) context%shape = shape
      if (present(estimate_shape)) context%estimate_shape = estimate_shape
      if (context%distribution /= dist_student) context%estimate_shape = .false.
      allocate(context%z(size(z,1),size(z,2)))
      context%z = z
      fit%distribution = context%distribution
      fit%time_varying = context%time_varying
      fit%shape = context%shape
      maxit = 800
      if (present(max_iterations)) maxit = max_iterations

      if (.not. context%time_varying .and. .not. context%estimate_shape) then
         allocate(fit%correlation(size(z,2),size(z,2)))
         fit%correlation = make_positive_definite(normalize_covariance(covariance_matrix(z)),1.0e-10_dp)
         fit%log_likelihood = static_copula_log_likelihood(z,fit%correlation, &
            context%distribution,context%shape,ok)
         fit%status = merge(0,2,ok)
         fit%iterations = 0
         k = size(z,2)*(size(z,2)-1)/2
      else
         ncoef = merge(context%p+context%q+context%g,0,context%time_varying)
         if (context%time_varying .and. ncoef == 0) then
            fit%status = 3
            return
         end if
         allocate(x0(ncoef+merge(1,0,context%estimate_shape)))
         if (ncoef > 0) then
            x0 = 0.0_dp
            if (context%p > 0) x0(1:context%p) = log((0.05_dp/real(context%p,dp))/0.03_dp)
            if (context%q > 0) x0(context%p+1:context%p+context%q) = &
               log((0.90_dp/real(context%q,dp))/0.03_dp)
            if (context%g > 0) x0(context%p+context%q+1:ncoef) = &
               log((0.01_dp/real(context%g,dp))/0.03_dp)
         end if
         if (context%estimate_shape) then
            fitted_shape = min(49.5_dp,max(2.1_dp,context%shape))
            remaining = (fitted_shape-2.01_dp)/(50.0_dp-fitted_shape)
            x0(size(x0)) = log(max(remaining,1.0e-8_dp))
         end if
         opt = nelder_mead(copula_objective,x0,context,step=0.16_dp,tolerance=1.0e-8_dp, &
            max_iterations=maxit)
         call decode_copula_parameters(opt%x,context,alpha,beta,gamma,fitted_shape)
         fit%shape = fitted_shape
         if (context%time_varying) then
            fit%dcc%spec = make_dcc_spec(alpha,beta,gamma,context%distribution,fitted_shape)
            allocate(fit%dcc%qbar(size(z,2),size(z,2)),fit%dcc%nbar(size(z,2),size(z,2)))
            allocate(fit%dcc%q(size(z,2),size(z,2),size(z,1)), &
               fit%dcc%r(size(z,2),size(z,2),size(z,1)),fit%dcc%loglikelihoods(size(z,1)), &
               fit%dcc%standardized_residuals(size(z,1),size(z,2)))
            fit%dcc%standardized_residuals = z
            call dcc_filter(z,fit%dcc%spec,fit%dcc%q,fit%dcc%r,fit%dcc%loglikelihoods, &
               fit%dcc%qbar,fit%dcc%nbar,ok)
            if (ok) then
               if (context%distribution == dist_gaussian) then
                  fit%log_likelihood = dynamic_gaussian_copula_log_likelihood(z,fit%dcc%spec,ok)
               else
                  fit%log_likelihood = dynamic_student_copula_log_likelihood(z,fit%dcc%spec, &
                     fitted_shape,ok)
               end if
            end if
            if (ok) then
               allocate(fit%correlation(size(z,2),size(z,2)))
               fit%correlation = fit%dcc%r(:,:,size(z,1))
               fit%dcc%log_likelihood = fit%log_likelihood
               fit%dcc%status = merge(0,opt%status,opt%status == 0)
               fit%dcc%iterations = opt%iterations
            end if
         else
            allocate(fit%correlation(size(z,2),size(z,2)))
            fit%correlation = make_positive_definite(normalize_covariance(covariance_matrix(z)),1.0e-10_dp)
            fit%log_likelihood = static_copula_log_likelihood(z,fit%correlation, &
               context%distribution,fitted_shape,ok)
         end if
         fit%iterations = opt%iterations
         fit%status = merge(0,2,ok)
         if (opt%status /= 0 .and. fit%status == 0) fit%status = opt%status
         k = size(opt%x)+size(z,2)*(size(z,2)-1)/2
      end if
      if (fit%status == 0 .or. fit%status == 1) then
         fit%aic = -2.0_dp*fit%log_likelihood+2.0_dp*real(k,dp)
         fit%bic = -2.0_dp*fit%log_likelihood+log(real(size(z,1),dp))*real(k,dp)
      end if
   end function fit_copula

   subroutine simulate_static_copula(nobs, correlation, distribution, scores, uniforms, shape, valid)
      integer, intent(in) :: nobs, distribution
      real(dp), intent(in) :: correlation(:,:)
      real(dp), intent(out) :: scores(nobs,size(correlation,1)), uniforms(nobs,size(correlation,1))
      real(dp), intent(in), optional :: shape
      logical, intent(out), optional :: valid
      real(dp) :: zero(size(correlation,1)), draw(size(correlation,1)), nu
      logical :: ok, step_ok
      integer :: t

      zero = 0.0_dp
      nu = 8.0_dp
      if (present(shape)) nu = shape
      ok = distribution == dist_gaussian .or. distribution == dist_student
      if (ok) then
         do t = 1, nobs
            if (distribution == dist_gaussian) then
               call random_multivariate_normal(zero,correlation,draw,step_ok)
            else
               call random_multivariate_student(zero,correlation,nu,draw,step_ok)
            end if
            if (step_ok) scores(t,:) = draw
            ok = ok .and. step_ok
         end do
         call score_to_uniform_transform(scores,distribution,uniforms,nu)
      else
         scores = 0.0_dp
         uniforms = 0.5_dp
      end if
      if (present(valid)) valid = ok
   end subroutine simulate_static_copula

   subroutine simulate_dynamic_copula(nobs, spec, qbar, scores, uniforms, q, r, burn)
      integer, intent(in) :: nobs
      type(dcc_spec), intent(in) :: spec
      real(dp), intent(in) :: qbar(:,:)
      real(dp), intent(out) :: scores(nobs,size(qbar,1)), uniforms(nobs,size(qbar,1))
      real(dp), intent(out) :: q(size(qbar,1),size(qbar,2),nobs)
      real(dp), intent(out) :: r(size(qbar,1),size(qbar,2),nobs)
      integer, intent(in), optional :: burn

      if (present(burn)) then
         call simulate_dcc(nobs,spec,qbar,scores,q,r,burn)
      else
         call simulate_dcc(nobs,spec,qbar,scores,q,r)
      end if
      call score_to_uniform_transform(scores,spec%distribution,uniforms,spec%shape)
   end subroutine simulate_dynamic_copula

   subroutine parametric_uniform_transform(x, distributions, shapes, u)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: distributions(:)
      real(dp), intent(in) :: shapes(:)
      real(dp), intent(out) :: u(size(x,1),size(x,2))
      integer :: i, j

      if (size(distributions) /= size(x,2) .or. size(shapes) /= size(x,2)) then
         u = 0.5_dp
         return
      end if
      do j = 1, size(x,2)
         do i = 1, size(x,1)
            select case (distributions(j))
            case (dist_gaussian)
               u(i,j) = normal_cdf(x(i,j))
            case (dist_student)
               u(i,j) = standardized_student_cdf(x(i,j),shapes(j))
            case (dist_laplace)
               u(i,j) = standardized_laplace_cdf(x(i,j))
            case default
               u(i,j) = 0.5_dp
            end select
            u(i,j) = clamp_probabilities(u(i,j))
         end do
      end do
   end subroutine parametric_uniform_transform

   subroutine copula_score_transform(u, distribution, z, shape)
      real(dp), intent(in) :: u(:,:)
      integer, intent(in) :: distribution
      real(dp), intent(out) :: z(size(u,1),size(u,2))
      real(dp), intent(in), optional :: shape
      real(dp) :: nu

      nu = 8.0_dp
      if (present(shape)) nu = shape
      select case (distribution)
      case (dist_gaussian)
         z = normal_quantile(clamp_probabilities(u))
      case (dist_student)
         z = standardized_student_quantile(clamp_probabilities(u),nu)
      case default
         z = 0.0_dp
      end select
   end subroutine copula_score_transform

   subroutine score_to_uniform_transform(z, distribution, u, shape)
      real(dp), intent(in) :: z(:,:)
      integer, intent(in) :: distribution
      real(dp), intent(out) :: u(size(z,1),size(z,2))
      real(dp), intent(in), optional :: shape
      real(dp) :: nu

      nu = 8.0_dp
      if (present(shape)) nu = shape
      select case (distribution)
      case (dist_gaussian)
         u = normal_cdf(z)
      case (dist_student)
         u = standardized_student_cdf(z,nu)
      case default
         u = 0.5_dp
      end select
      u = clamp_probabilities(u)
   end subroutine score_to_uniform_transform

   function copula_objective(x, generic_context) result(value)
      real(dp), intent(in) :: x(:)
      class(*), intent(in) :: generic_context
      real(dp) :: value, shape, ll
      real(dp), allocatable :: alpha(:), beta(:), gamma(:)
      type(dcc_spec) :: spec
      real(dp), allocatable :: rstatic(:,:)
      logical :: ok

      select type (context => generic_context)
      type is (copula_context)
         call decode_copula_parameters(x,context,alpha,beta,gamma,shape)
         if (context%time_varying) then
            spec = make_dcc_spec(alpha,beta,gamma,context%distribution,shape)
            if (context%distribution == dist_gaussian) then
               ll = dynamic_gaussian_copula_log_likelihood(context%z,spec,ok)
            else
               ll = dynamic_student_copula_log_likelihood(context%z,spec,shape,ok)
            end if
         else
            allocate(rstatic(size(context%z,2),size(context%z,2)))
            rstatic = make_positive_definite(normalize_covariance(covariance_matrix(context%z)),1.0e-10_dp)
            ll = static_copula_log_likelihood(context%z,rstatic,context%distribution,shape,ok)
         end if
         if (ok) then
            value = -ll
         else
            value = huge(1.0_dp)/100.0_dp
         end if
      class default
         value = huge(1.0_dp)/100.0_dp
      end select
   end function copula_objective

   subroutine decode_copula_parameters(x, context, alpha, beta, gamma, shape)
      real(dp), intent(in) :: x(:)
      type(copula_context), intent(in) :: context
      real(dp), allocatable, intent(out) :: alpha(:), beta(:), gamma(:)
      real(dp), intent(out) :: shape
      real(dp), allocatable :: ex(:), weights(:)
      real(dp) :: denominator, logistic
      integer :: ncoef, offset

      ncoef = merge(context%p+context%q+context%g,0,context%time_varying)
      allocate(alpha(context%p),beta(context%q),gamma(context%g))
      if (ncoef > 0) then
         allocate(ex(ncoef),weights(ncoef))
         ex = exp(max(-30.0_dp,min(30.0_dp,x(1:ncoef))))
         denominator = 1.0_dp+sum(ex)
         weights = 0.999_dp*ex/denominator
         if (context%p > 0) alpha = weights(1:context%p)
         offset = context%p
         if (context%q > 0) beta = weights(offset+1:offset+context%q)
         offset = offset+context%q
         if (context%g > 0) gamma = weights(offset+1:offset+context%g)
      else
         alpha = 0.0_dp
         beta = 0.0_dp
         gamma = 0.0_dp
      end if
      shape = context%shape
      if (context%estimate_shape) then
         logistic = 1.0_dp/(1.0_dp+exp(-max(-30.0_dp,min(30.0_dp,x(size(x))))))
         shape = 2.01_dp+47.99_dp*logistic
      end if
   end subroutine decode_copula_parameters

end module rmgarch_copula
