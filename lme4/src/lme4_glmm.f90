module lme4_glmm
   use lme4_kinds, only : dp, pi
   use lme4_types, only : random_term_t, covariance_block_t, glmm_control_t, &
      glmm_result_t, family_binomial, family_poisson, family_gamma, &
      family_inverse_gaussian, family_negative_binomial, &
      covariance_unstructured, covariance_diagonal, covariance_compound_symmetry, covariance_ar1
   use lme4_covariance, only : validate_terms, total_theta, build_random_design, &
      build_covariance_from_eta, eta_to_theta
   use lme4_linalg, only : cholesky_lower, chol_solve, invert_spd, logdet_from_chol
   use minqa_module, only : bobyqa, minqa_control_t, minqa_result_t
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private
   public :: fit_glmm, fit_glmer_nb, predict_glmm

   real(dp), allocatable :: active_y(:), active_x(:,:), active_z(:,:)
   real(dp), allocatable :: active_weights(:), active_offset(:)
   type(random_term_t), allocatable :: active_terms(:)
   type(glmm_control_t) :: active_control
   integer :: active_family = 0
   real(dp) :: active_dispersion = 1.0_dp

contains

   subroutine fit_glmm(y, x, terms, family, result, weights, offset, control, dispersion)
      real(dp), intent(in) :: y(:), x(:,:)
      type(random_term_t), intent(inout) :: terms(:)
      integer, intent(in) :: family
      type(glmm_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), offset(:), dispersion
      type(glmm_control_t), intent(in), optional :: control

      type(glmm_control_t) :: ctrl
      type(minqa_control_t) :: mctrl
      type(minqa_result_t) :: mresult
      real(dp), allocatable :: eta_cov(:), lower(:), upper(:), z(:,:), w(:), off(:)
      real(dp), allocatable :: beta(:), u(:), linear_predictor(:), mu(:), residuals(:)
      real(dp), allocatable :: vcov_beta(:,:)
      type(covariance_block_t), allocatable :: blocks(:)
      integer, allocatable :: offsets(:)
      real(dp) :: criterion, disp
      integer :: n, p, nt, info, evals, pirls_iterations
      logical :: ok
      character(len=:), allocatable :: message

      ctrl = glmm_control_t()
      if (present(control)) ctrl = control
      disp = 1.0_dp
      if (present(dispersion)) disp = dispersion
      call initialize_result(result, family)
      result%dispersion = disp
      n = size(y)
      p = size(x,2)
      if (size(x,1) /= n .or. n <= p .or. p < 1) then
         call fail_result(result, 1, &
            'require size(x,1)=size(y), at least one fixed effect, and n>p')
         return
      end if
      if (.not. supported_family(family)) then
         call fail_result(result, 1, &
            'supported GLMM families are binomial, poisson, gamma, inverse Gaussian, and negative binomial')
         return
      end if
      if ((family == family_gamma .or. family == family_inverse_gaussian .or. &
           family == family_negative_binomial) .and. disp <= 0.0_dp) then
         call fail_result(result, 1, 'dispersion or negative-binomial size must be positive')
         return
      end if
      call validate_response(family, y, ok, message)
      if (.not. ok) then
         call fail_result(result, 1, message)
         return
      end if
      call validate_terms(terms, n, ok, message)
      if (.not. ok) then
         call fail_result(result, 1, message)
         return
      end if
      allocate(w(n), off(n))
      w = 1.0_dp
      off = 0.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights <= 0.0_dp)) then
            call fail_result(result, 1, 'weights must be positive and match the response length')
            return
         end if
         w = weights
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            call fail_result(result, 1, 'offset must match the response length')
            return
         end if
         off = offset
      end if
      call build_random_design(terms, z, offsets)
      nt = total_theta(terms)
      allocate(eta_cov(nt), lower(nt), upper(nt))
      call initialize_parameters(terms, eta_cov, lower, upper, ctrl)

      active_y = y
      active_x = x
      active_z = z
      active_weights = w
      active_offset = off
      active_terms = terms
      active_family = family
      active_dispersion = disp
      active_control = ctrl

      if (nt == 1) then
         call golden_minimize(lower(1), upper(1), ctrl%tolerance, ctrl%maxfun, &
            eta_cov(1), criterion, evals)
         info = 0
      else
         mctrl%maxfun = ctrl%maxfun
         mctrl%rhoend = max(1.0e-9_dp, ctrl%tolerance)
         mctrl%rhobeg = 0.25_dp
         mctrl%npt = min((nt+1)*(nt+2)/2, max(nt+2, 2*nt+1))
         call bobyqa(glmm_objective, eta_cov, mresult, lower, upper, mctrl)
         criterion = mresult%fval
         evals = mresult%evaluations
         info = mresult%status
      end if

      call evaluate_glmm(eta_cov, criterion, beta, u, linear_predictor, mu, &
         residuals, vcov_beta, blocks, pirls_iterations, info)
      if (info /= 0) then
         call fail_result(result, 2, 'GLMM optimization or final PIRLS iteration failed')
         call clear_active()
         return
      end if
      result%beta = beta
      result%u = u
      call eta_to_theta(terms, eta_cov, result%theta)
      result%linear_predictor = linear_predictor
      result%fitted = mu
      result%residuals = residuals
      result%vcov_beta = vcov_beta
      result%varcorr = blocks
      result%term_offsets = offsets
      result%deviance = criterion
      result%log_likelihood = -0.5_dp*criterion
      result%aic = criterion + 2.0_dp*real(p+nt,dp)
      result%bic = criterion + log(real(n,dp))*real(p+nt,dp)
      result%dispersion = disp
      result%quadrature_order = 1
      result%evaluations = evals
      result%pirls_iterations = pirls_iterations
      result%status = 0
      result%converged = .true.
      result%message = 'converged'
      call clear_active()
   end subroutine fit_glmm

   subroutine fit_glmer_nb(y, x, terms, result, weights, offset, control)
      real(dp), intent(in) :: y(:), x(:,:)
      type(random_term_t), intent(inout) :: terms(:)
      type(glmm_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), offset(:)
      type(glmm_control_t), intent(in), optional :: control

      type(glmm_control_t) :: ctrl
      type(glmm_result_t) :: trial
      real(dp), allocatable :: w(:), off(:)
      real(dp), parameter :: gr = 0.6180339887498948482_dp
      real(dp) :: left, right, c, d, fc, fd, logtheta
      integer :: evaluations

      ctrl = glmm_control_t()
      if (present(control)) ctrl = control
      allocate(w(size(y)), off(size(y)))
      w = 1.0_dp
      off = 0.0_dp
      if (present(weights)) then
         if (size(weights) /= size(y)) then
            call initialize_result(result, family_negative_binomial)
            call fail_result(result, 1, 'weights must match the response length')
            return
         end if
         w = weights
      end if
      if (present(offset)) then
         if (size(offset) /= size(y)) then
            call initialize_result(result, family_negative_binomial)
            call fail_result(result, 1, 'offset must match the response length')
            return
         end if
         off = offset
      end if

      left = ctrl%lower_log_dispersion
      right = ctrl%upper_log_dispersion
      c = right-gr*(right-left)
      d = left+gr*(right-left)
      call profile_value(c, fc)
      call profile_value(d, fd)
      evaluations = 2
      do while (abs(right-left) > ctrl%tolerance*(1.0_dp+abs(left)+abs(right)) .and. &
                evaluations < ctrl%max_profile)
         if (fc < fd) then
            right = d
            d = c
            fd = fc
            c = right-gr*(right-left)
            call profile_value(c, fc)
         else
            left = c
            c = d
            fc = fd
            d = left+gr*(right-left)
            call profile_value(d, fd)
         end if
         evaluations = evaluations+1
      end do
      if (fc < fd) then
         logtheta = c
      else
         logtheta = d
      end if
      call fit_glmm(y, x, terms, family_negative_binomial, result, weights=w, &
         offset=off, control=ctrl, dispersion=exp(logtheta))
      if (result%converged) then
         result%aic = result%aic + 2.0_dp
         result%bic = result%bic + log(real(size(y),dp))
         result%evaluations = result%evaluations + evaluations
         result%message = 'converged with profiled negative-binomial size'
      end if

   contains

      subroutine profile_value(log_size, value)
         real(dp), intent(in) :: log_size
         real(dp), intent(out) :: value
         call fit_glmm(y, x, terms, family_negative_binomial, trial, weights=w, &
            offset=off, control=ctrl, dispersion=exp(log_size))
         if (trial%converged .and. ieee_is_finite(trial%deviance)) then
            value = trial%deviance
         else
            value = huge(1.0_dp)/100.0_dp
         end if
      end subroutine profile_value

   end subroutine fit_glmer_nb

   real(dp) function glmm_objective(eta_cov) result(value)
      real(dp), intent(in) :: eta_cov(:)
      integer :: info, iterations
      call evaluate_glmm(eta_cov, value, pirls_iterations=iterations, info=info)
      if (info /= 0 .or. .not. ieee_is_finite(value)) value = huge(1.0_dp)/100.0_dp
   end function glmm_objective

   subroutine evaluate_glmm(eta_cov, criterion, beta, u, linear_predictor, mu, &
      residuals, vcov_beta, blocks, pirls_iterations, info)
      real(dp), intent(in) :: eta_cov(:)
      real(dp), intent(out) :: criterion
      real(dp), allocatable, intent(out), optional :: beta(:), u(:)
      real(dp), allocatable, intent(out), optional :: linear_predictor(:), mu(:), residuals(:)
      real(dp), allocatable, intent(out), optional :: vcov_beta(:,:)
      type(covariance_block_t), allocatable, intent(out), optional :: blocks(:)
      integer, intent(out) :: pirls_iterations
      integer, intent(out) :: info

      real(dp), allocatable :: g(:,:), ginv(:,:), h(:,:), lh(:,:)
      real(dp), allocatable :: b(:), raneff(:), eta(:), mean(:), work_w(:), work_z(:)
      real(dp), allocatable :: a(:,:), la(:,:), rhs(:), coef(:), oldcoef(:), target(:)
      real(dp), allocatable :: fullinv(:,:), covb(:,:)
      type(covariance_block_t), allocatable :: relblocks(:)
      real(dp) :: logdetg, logdeth, penalty, neg2ll, change
      integer :: n, p, nr, iter, covinfo

      criterion = huge(1.0_dp)/100.0_dp
      info = 0
      n = size(active_y)
      p = size(active_x,2)
      nr = size(active_z,2)
      call build_covariance_from_eta(active_terms, eta_cov, g, relblocks, covinfo)
      if (covinfo /= 0) then
         info = 1
         pirls_iterations = 0
         return
      end if
      call invert_spd(g, ginv, info, logdetg)
      if (info /= 0) then
         pirls_iterations = 0
         return
      end if
      allocate(b(p), raneff(nr), eta(n), mean(n), work_w(n), work_z(n))
      allocate(coef(p+nr), oldcoef(p+nr), target(n))
      b = 0.0_dp
      raneff = 0.0_dp
      if (active_family == family_binomial) then
         b(1) = log(max(1.0e-6_dp,min(1.0_dp-1.0e-6_dp, &
            sum(active_weights*active_y)/sum(active_weights)))/ &
            max(1.0e-6_dp,1.0_dp-sum(active_weights*active_y)/sum(active_weights)))
      else
         b(1) = log(max(1.0e-6_dp,sum(active_weights*active_y)/sum(active_weights)))
      end if
      coef(1:p) = b
      coef(p+1:p+nr) = raneff

      do iter = 1, active_control%max_pirls
         oldcoef = coef
         eta = active_offset + matmul(active_x,coef(1:p)) + &
            matmul(active_z,coef(p+1:p+nr))
         call family_working(active_family, active_y, eta, active_weights, &
            active_dispersion, mean, work_w, work_z)
         target = work_z-active_offset
         allocate(a(p+nr,p+nr), rhs(p+nr))
         call weighted_cross_system(active_x, active_z, work_w, target, ginv, a, rhs)
         call cholesky_lower(a, la, info, jitter=1.0e-10_dp)
         if (info /= 0) then
            pirls_iterations = iter
            return
         end if
         call chol_solve(la, rhs, coef)
         change = maxval(abs(coef-oldcoef)/(1.0_dp+abs(oldcoef)))
         deallocate(a, rhs, la)
         if (change < active_control%pirls_tolerance) exit
      end do
      pirls_iterations = iter
      if (iter > active_control%max_pirls) then
         info = 3
         return
      end if
      b = coef(1:p)
      raneff = coef(p+1:p+nr)
      eta = active_offset+matmul(active_x,b)+matmul(active_z,raneff)
      call family_working(active_family, active_y, eta, active_weights, &
         active_dispersion, mean, work_w, work_z)
      h = matmul(transpose(active_z),spread(work_w,2,nr)*active_z)+ginv
      call cholesky_lower(h, lh, info)
      if (info /= 0) return
      logdeth = logdet_from_chol(lh)
      penalty = dot_product(raneff,matmul(ginv,raneff))
      neg2ll = conditional_deviance(active_family,active_y,mean,active_weights, &
         active_dispersion)
      criterion = neg2ll+penalty+logdetg+logdeth

      if (present(beta)) beta = b
      if (present(u)) u = raneff
      if (present(linear_predictor)) linear_predictor = eta
      if (present(mu)) mu = mean
      if (present(residuals)) residuals = active_y-mean
      if (present(vcov_beta)) then
         allocate(a(p+nr,p+nr),rhs(p+nr))
         target = work_z-active_offset
         call weighted_cross_system(active_x,active_z,work_w,target,ginv,a,rhs)
         call invert_spd(a,fullinv,info)
         if (info /= 0) return
         allocate(covb(p,p))
         covb = fullinv(1:p,1:p)
         vcov_beta = covb
      end if
      if (present(blocks)) blocks = relblocks
   end subroutine evaluate_glmm

   subroutine family_working(family, y, eta, prior_weights, dispersion, mean, work_w, work_z)
      integer, intent(in) :: family
      real(dp), intent(in) :: y(:), eta(:), prior_weights(:), dispersion
      real(dp), intent(out) :: mean(:), work_w(:), work_z(:)
      integer :: i
      real(dp) :: derivative
      do i = 1, size(y)
         select case (family)
         case (family_binomial)
            mean(i) = logistic(eta(i))
            mean(i) = min(1.0_dp-1.0e-10_dp,max(1.0e-10_dp,mean(i)))
            derivative = mean(i)*(1.0_dp-mean(i))
            work_w(i) = max(1.0e-10_dp,prior_weights(i)*derivative)
            work_z(i) = eta(i)+(y(i)-mean(i))/derivative
         case (family_poisson)
            mean(i) = exp(min(30.0_dp,max(-30.0_dp,eta(i))))
            derivative = max(1.0e-10_dp,mean(i))
            work_w(i) = max(1.0e-10_dp,prior_weights(i)*mean(i))
            work_z(i) = eta(i)+(y(i)-mean(i))/derivative
         case (family_gamma)
            mean(i) = exp(min(30.0_dp,max(-30.0_dp,eta(i))))
            derivative = max(1.0e-10_dp,mean(i))
            work_w(i) = max(1.0e-10_dp,prior_weights(i)/dispersion)
            work_z(i) = eta(i)+(y(i)-mean(i))/derivative
         case (family_inverse_gaussian)
            mean(i) = exp(min(30.0_dp,max(-30.0_dp,eta(i))))
            derivative = max(1.0e-10_dp,mean(i))
            work_w(i) = max(1.0e-10_dp,prior_weights(i)/(dispersion*mean(i)))
            work_z(i) = eta(i)+(y(i)-mean(i))/derivative
         case (family_negative_binomial)
            mean(i) = exp(min(30.0_dp,max(-30.0_dp,eta(i))))
            derivative = max(1.0e-10_dp,mean(i))
            work_w(i) = max(1.0e-10_dp,prior_weights(i)*mean(i)/ &
               (1.0_dp+mean(i)/dispersion))
            work_z(i) = eta(i)+(y(i)-mean(i))/derivative
         end select
      end do
   end subroutine family_working

   subroutine weighted_cross_system(x, z, w, target, ginv, a, rhs)
      real(dp), intent(in) :: x(:,:), z(:,:), w(:), target(:), ginv(:,:)
      real(dp), intent(out) :: a(:,:), rhs(:)
      real(dp), allocatable :: d(:,:)
      integer :: p, nr, i
      p = size(x,2)
      nr = size(z,2)
      allocate(d(size(x,1),p+nr))
      d(:,1:p) = x
      d(:,p+1:p+nr) = z
      a = matmul(transpose(d),spread(w,2,p+nr)*d)
      a(p+1:p+nr,p+1:p+nr) = a(p+1:p+nr,p+1:p+nr)+ginv
      do i = 1, p
         a(i,i) = a(i,i)+1.0e-12_dp
      end do
      rhs = matmul(transpose(d),w*target)
   end subroutine weighted_cross_system

   real(dp) function conditional_deviance(family,y,mean,weights,dispersion) result(value)
      integer, intent(in) :: family
      real(dp), intent(in) :: y(:), mean(:), weights(:), dispersion
      integer :: i
      real(dp) :: shape, lambda, loglik
      value = 0.0_dp
      select case (family)
      case (family_binomial)
         do i = 1, size(y)
            value = value-2.0_dp*weights(i)*(y(i)*log(mean(i)) + &
               (1.0_dp-y(i))*log(1.0_dp-mean(i)))
         end do
      case (family_poisson)
         do i = 1, size(y)
            value = value-2.0_dp*weights(i)*(y(i)*log(mean(i))-mean(i)- &
               log_gamma(y(i)+1.0_dp))
         end do
      case (family_gamma)
         shape = 1.0_dp/dispersion
         do i = 1, size(y)
            loglik = (shape-1.0_dp)*log(y(i))-y(i)/(dispersion*mean(i)) - &
               shape*log(dispersion*mean(i))-log_gamma(shape)
            value = value-2.0_dp*weights(i)*loglik
         end do
      case (family_inverse_gaussian)
         lambda = 1.0_dp/dispersion
         do i = 1, size(y)
            loglik = 0.5_dp*(log(lambda)-log(2.0_dp*pi)-3.0_dp*log(y(i))) - &
               lambda*(y(i)-mean(i))**2/(2.0_dp*mean(i)**2*y(i))
            value = value-2.0_dp*weights(i)*loglik
         end do
      case (family_negative_binomial)
         do i = 1, size(y)
            loglik = log_gamma(y(i)+dispersion)-log_gamma(dispersion) - &
               log_gamma(y(i)+1.0_dp)+dispersion*log(dispersion/(dispersion+mean(i))) + &
               y(i)*log(mean(i)/(dispersion+mean(i)))
            value = value-2.0_dp*weights(i)*loglik
         end do
      end select
   end function conditional_deviance

   pure logical function supported_family(family) result(ok)
      integer, intent(in) :: family
      ok = family == family_binomial .or. family == family_poisson .or. &
         family == family_gamma .or. family == family_inverse_gaussian .or. &
         family == family_negative_binomial
   end function supported_family

   subroutine validate_response(family, y, ok, message)
      integer, intent(in) :: family
      real(dp), intent(in) :: y(:)
      logical, intent(out) :: ok
      character(len=:), allocatable, intent(out) :: message
      ok = .false.
      select case (family)
      case (family_binomial)
         if (any(y < 0.0_dp) .or. any(y > 1.0_dp)) then
            message = 'binomial responses must lie in [0,1]'
            return
         end if
      case (family_poisson, family_negative_binomial)
         if (any(y < 0.0_dp)) then
            message = 'count responses must be nonnegative'
            return
         end if
      case (family_gamma, family_inverse_gaussian)
         if (any(y <= 0.0_dp)) then
            message = 'gamma and inverse-Gaussian responses must be positive'
            return
         end if
      end select
      ok = .true.
      message = 'ok'
   end subroutine validate_response

   pure real(dp) function logistic(x) result(value)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         value = 1.0_dp/(1.0_dp+exp(-min(700.0_dp,x)))
      else
         value = exp(max(-700.0_dp,x))/(1.0_dp+exp(max(-700.0_dp,x)))
      end if
   end function logistic

   subroutine initialize_parameters(terms, eta_cov, lower, upper, ctrl)
      type(random_term_t), intent(in) :: terms(:)
      real(dp), intent(out) :: eta_cov(:), lower(:), upper(:)
      type(glmm_control_t), intent(in) :: ctrl
      integer :: k, q, i, j, idx

      idx = 0
      do k = 1, size(terms)
         q = terms(k)%n_coefficients()
         select case (terms(k)%covariance_structure)
         case (covariance_unstructured)
            do j = 1, q
               do i = j, q
                  idx = idx+1
                  if (i == j) then
                     eta_cov(idx) = log(0.5_dp)
                     lower(idx) = ctrl%lower_log_sd
                     upper(idx) = ctrl%upper_log_sd
                  else
                     eta_cov(idx) = 0.0_dp
                     lower(idx) = ctrl%lower_offdiag
                     upper(idx) = ctrl%upper_offdiag
                  end if
               end do
            end do
         case (covariance_diagonal)
            do i = 1, q
               idx = idx+1
               eta_cov(idx) = log(0.5_dp)
               lower(idx) = ctrl%lower_log_sd
               upper(idx) = ctrl%upper_log_sd
            end do
         case (covariance_compound_symmetry, covariance_ar1)
            idx = idx+1
            eta_cov(idx) = log(0.5_dp)
            lower(idx) = ctrl%lower_log_sd
            upper(idx) = ctrl%upper_log_sd
            idx = idx+1
            eta_cov(idx) = 0.0_dp
            lower(idx) = ctrl%lower_offdiag
            upper(idx) = ctrl%upper_offdiag
         end select
      end do
   end subroutine initialize_parameters

   subroutine golden_minimize(a, b, tol, maxfun, xmin, fmin, evaluations)
      real(dp), intent(in) :: a, b, tol
      integer, intent(in) :: maxfun
      real(dp), intent(out) :: xmin, fmin
      integer, intent(out) :: evaluations
      real(dp), parameter :: gr = 0.6180339887498948482_dp
      real(dp) :: left, right, c, d, fc, fd
      real(dp) :: xv(1)
      left = a
      right = b
      c = right-gr*(right-left)
      d = left+gr*(right-left)
      xv(1) = c
      fc = glmm_objective(xv)
      xv(1) = d
      fd = glmm_objective(xv)
      evaluations = 2
      do while (abs(right-left) > tol*(1.0_dp+abs(left)+abs(right)) .and. &
                evaluations < maxfun)
         if (fc < fd) then
            right = d
            d = c
            fd = fc
            c = right-gr*(right-left)
            xv(1) = c
            fc = glmm_objective(xv)
         else
            left = c
            c = d
            fc = fd
            d = left+gr*(right-left)
            xv(1) = d
            fd = glmm_objective(xv)
         end if
         evaluations = evaluations+1
      end do
      if (fc < fd) then
         xmin = c
         fmin = fc
      else
         xmin = d
         fmin = fd
      end if
   end subroutine golden_minimize

   subroutine predict_glmm(result, x_new, z_new, prediction, response_scale, offset)
      type(glmm_result_t), intent(in) :: result
      real(dp), intent(in) :: x_new(:,:), z_new(:,:)
      real(dp), allocatable, intent(out) :: prediction(:)
      logical, intent(in), optional :: response_scale
      real(dp), intent(in), optional :: offset(:)
      logical :: response
      integer :: i
      response = .true.
      if (present(response_scale)) response = response_scale
      prediction = matmul(x_new,result%beta)+matmul(z_new,result%u)
      if (present(offset)) prediction = prediction+offset
      if (response) then
         select case (result%family)
         case (family_binomial)
            do i = 1, size(prediction)
               prediction(i) = logistic(prediction(i))
            end do
         case (family_poisson, family_gamma, family_inverse_gaussian, &
               family_negative_binomial)
            prediction = exp(min(30.0_dp,max(-30.0_dp,prediction)))
         end select
      end if
   end subroutine predict_glmm

   subroutine initialize_result(result, family)
      type(glmm_result_t), intent(out) :: result
      integer, intent(in) :: family
      result%family = family
      result%message = 'not fitted'
   end subroutine initialize_result

   subroutine fail_result(result, status, message)
      type(glmm_result_t), intent(inout) :: result
      integer, intent(in) :: status
      character(len=*), intent(in) :: message
      result%status = status
      result%converged = .false.
      result%message = message
   end subroutine fail_result

   subroutine clear_active()
      if (allocated(active_y)) deallocate(active_y)
      if (allocated(active_x)) deallocate(active_x)
      if (allocated(active_z)) deallocate(active_z)
      if (allocated(active_weights)) deallocate(active_weights)
      if (allocated(active_offset)) deallocate(active_offset)
      if (allocated(active_terms)) deallocate(active_terms)
      active_dispersion = 1.0_dp
   end subroutine clear_active

end module lme4_glmm
