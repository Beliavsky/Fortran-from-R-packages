module lme4_aghq
   use lme4_kinds, only : dp, pi
   use lme4_types, only : random_term_t, covariance_block_t, glmm_control_t, &
      glmm_result_t, gh_rule_t, family_binomial, family_poisson, &
      family_negative_binomial
   use lme4_covariance, only : validate_terms
   use lme4_quadrature, only : gh_rule
   use lme4_linalg, only : invert_spd
   use lme4_glmm, only : fit_glmm
   use minqa_module, only : bobyqa, minqa_control_t, minqa_result_t
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private
   public :: fit_glmm_aghq

   real(dp), allocatable :: aghq_y(:), aghq_x(:,:), aghq_z(:)
   real(dp), allocatable :: aghq_weights(:), aghq_offset(:)
   integer, allocatable :: aghq_group(:)
   type(gh_rule_t) :: aghq_rule
   integer :: aghq_family = 0
   integer :: aghq_levels = 0
   real(dp) :: aghq_dispersion = 1.0_dp

contains

   subroutine fit_glmm_aghq(y, x, terms, family, result, order, weights, offset, &
      dispersion, control)
      real(dp), intent(in) :: y(:), x(:,:)
      type(random_term_t), intent(inout) :: terms(:)
      integer, intent(in) :: family, order
      type(glmm_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), offset(:), dispersion
      type(glmm_control_t), intent(in), optional :: control

      type(glmm_control_t) :: ctrl
      type(glmm_result_t) :: initial
      type(minqa_control_t) :: mctrl
      type(minqa_result_t) :: mresult
      real(dp), allocatable :: w(:), off(:), par(:), lower(:), upper(:)
      real(dp), allocatable :: modes(:), eta(:), mu(:), hess(:,:), hinv(:,:)
      real(dp) :: disp, criterion, sigma
      integer :: n, p, npar, info, i, hinfo
      logical :: ok
      character(len=:), allocatable :: message

      ctrl = glmm_control_t()
      if (present(control)) ctrl = control
      disp = 1.0_dp
      if (present(dispersion)) disp = dispersion
      call initialize_result(result, family, disp, order)
      n = size(y)
      p = size(x,2)
      if (order < 1) then
         call fail_result(result, 1, 'adaptive quadrature order must be positive')
         return
      end if
      if (size(x,1) /= n .or. p < 1 .or. n <= p) then
         call fail_result(result, 1, 'invalid fixed-effect design matrix')
         return
      end if
      if (family /= family_binomial .and. family /= family_poisson .and. &
          family /= family_negative_binomial) then
         call fail_result(result, 1, &
            'AGHQ supports binomial, poisson, and negative-binomial families')
         return
      end if
      if (family == family_binomial .and. (any(y < 0.0_dp) .or. any(y > 1.0_dp))) then
         call fail_result(result, 1, 'binomial responses must lie in [0,1]')
         return
      end if
      if ((family == family_poisson .or. family == family_negative_binomial) .and. &
          any(y < 0.0_dp)) then
         call fail_result(result, 1, 'count responses must be nonnegative')
         return
      end if
      if (family == family_negative_binomial .and. disp <= 0.0_dp) then
         call fail_result(result, 1, 'negative-binomial size must be positive')
         return
      end if
      call validate_terms(terms, n, ok, message)
      if (.not. ok) then
         call fail_result(result, 1, message)
         return
      end if
      if (size(terms) /= 1 .or. terms(1)%n_coefficients() /= 1) then
         call fail_result(result, 1, &
            'AGHQ currently requires one grouped scalar random-effect term')
         return
      end if

      allocate(w(n), off(n))
      w = 1.0_dp
      off = 0.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights <= 0.0_dp)) then
            call fail_result(result, 1, 'weights must be positive and match y')
            return
         end if
         w = weights
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            call fail_result(result, 1, 'offset must match y')
            return
         end if
         off = offset
      end if
      call gh_rule(order, aghq_rule, info)
      if (info /= 0) then
         call fail_result(result, 1, 'failed to construct Gauss-Hermite rule')
         return
      end if

      call fit_glmm(y, x, terms, family, initial, weights=w, offset=off, &
         control=ctrl, dispersion=disp)
      npar = p+1
      allocate(par(npar), lower(npar), upper(npar))
      if (initial%converged) then
         par(1:p) = initial%beta
         par(npar) = log(max(1.0e-8_dp,initial%theta(1)))
      else
         par = 0.0_dp
         if (family == family_binomial) then
            par(1) = log(max(1.0e-6_dp,min(1.0_dp-1.0e-6_dp,sum(w*y)/sum(w)))/ &
               max(1.0e-6_dp,1.0_dp-sum(w*y)/sum(w)))
         else
            par(1) = log(max(1.0e-6_dp,sum(w*y)/sum(w)))
         end if
         par(npar) = log(0.5_dp)
      end if
      lower(1:p) = ctrl%lower_beta
      upper(1:p) = ctrl%upper_beta
      lower(npar) = ctrl%lower_log_sd
      upper(npar) = ctrl%upper_log_sd
      par = min(upper,max(lower,par))

      aghq_y = y
      aghq_x = x
      aghq_z = terms(1)%z(:,1)
      aghq_group = terms(1)%group
      aghq_weights = w
      aghq_offset = off
      aghq_family = family
      aghq_dispersion = disp
      aghq_levels = terms(1)%n_levels

      mctrl%maxfun = ctrl%maxfun
      mctrl%rhoend = max(1.0e-8_dp,ctrl%tolerance)
      mctrl%rhobeg = min(0.25_dp,0.2_dp*minval(upper-lower))
      mctrl%npt = min((npar+1)*(npar+2)/2,max(npar+2,2*npar+1))
      call bobyqa(aghq_objective, par, mresult, lower, upper, mctrl)
      criterion = mresult%fval
      if (.not. ieee_is_finite(criterion)) then
         call fail_result(result, 2, 'AGHQ optimization failed')
         call clear_active()
         return
      end if
      sigma = exp(par(npar))
      allocate(modes(aghq_levels))
      call evaluate_modes(par(1:p),sigma,modes,criterion,info)
      if (info /= 0) then
         call fail_result(result, 2, 'failed to evaluate adaptive quadrature modes')
         call clear_active()
         return
      end if

      result%beta = par(1:p)
      result%u = modes
      allocate(result%theta(1))
      result%theta(1) = sigma
      allocate(result%varcorr(1))
      result%varcorr(1)%name = terms(1)%name
      result%varcorr(1)%n_levels = terms(1)%n_levels
      allocate(result%varcorr(1)%covariance(1,1),result%varcorr(1)%sdcor(1,1))
      result%varcorr(1)%covariance(1,1) = sigma*sigma
      result%varcorr(1)%sdcor(1,1) = sigma
      allocate(result%term_offsets(2))
      result%term_offsets = [1,terms(1)%n_levels+1]
      eta = off+matmul(x,result%beta)
      do i = 1, n
         eta(i) = eta(i)+terms(1)%z(i,1)*modes(terms(1)%group(i))
      end do
      allocate(mu(n))
      select case (family)
      case (family_binomial)
         do i = 1, n
            mu(i) = logistic(eta(i))
         end do
      case default
         mu = exp(min(30.0_dp,max(-30.0_dp,eta)))
      end select
      result%linear_predictor = eta
      result%fitted = mu
      result%residuals = y-mu
      result%deviance = criterion
      result%log_likelihood = -0.5_dp*criterion
      result%aic = criterion+2.0_dp*real(npar,dp)
      result%bic = criterion+log(real(n,dp))*real(npar,dp)
      result%dispersion = disp
      result%quadrature_order = order
      result%evaluations = mresult%evaluations
      result%pirls_iterations = 0
      result%status = 0
      result%converged = .true.
      result%message = 'converged with adaptive Gauss-Hermite quadrature'

      call numerical_hessian(par,hess)
      call invert_spd(hess,hinv,hinfo)
      allocate(result%vcov_beta(p,p))
      if (hinfo == 0) then
         result%vcov_beta = 2.0_dp*hinv(1:p,1:p)
      else
         result%vcov_beta = 0.0_dp
      end if
      call clear_active()
   end subroutine fit_glmm_aghq

   real(dp) function aghq_objective(par) result(value)
      real(dp), intent(in) :: par(:)
      real(dp), allocatable :: modes(:)
      real(dp) :: sigma
      integer :: info
      sigma = exp(par(size(par)))
      allocate(modes(aghq_levels))
      call evaluate_modes(par(1:size(par)-1),sigma,modes,value,info)
      if (info /= 0 .or. .not. ieee_is_finite(value)) value = huge(1.0_dp)/100.0_dp
   end function aghq_objective

   subroutine evaluate_modes(beta,sigma,modes,criterion,info)
      real(dp), intent(in) :: beta(:), sigma
      real(dp), intent(out) :: modes(:), criterion
      integer, intent(out) :: info
      real(dp) :: log_integral
      integer :: g
      criterion = 0.0_dp
      info = 0
      do g = 1, aghq_levels
         call cluster_mode_integral(beta,sigma,g,modes(g),log_integral,info)
         if (info /= 0) return
         criterion = criterion-2.0_dp*log_integral
      end do
   end subroutine evaluate_modes

   subroutine cluster_mode_integral(beta,sigma,level,mode,log_integral,info)
      real(dp), intent(in) :: beta(:), sigma
      integer, intent(in) :: level
      real(dp), intent(out) :: mode, log_integral
      integer, intent(out) :: info
      real(dp), allocatable :: log_terms(:)
      real(dp) :: h, grad, hess, trial_h, trial_grad, trial_hess
      real(dp) :: step, trial, scale, maximum
      integer :: iter, j, shrink

      mode = 0.0_dp
      info = 0
      call cluster_derivatives(beta,sigma,level,mode,h,grad,hess)
      do iter = 1, 60
         if (hess >= -1.0e-12_dp) then
            info = 1
            return
         end if
         step = -grad/hess
         step = sign(min(abs(step),max(1.0_dp,2.0_dp*sigma)),step)
         trial = mode+step
         do shrink = 1, 20
            call cluster_derivatives(beta,sigma,level,trial,trial_h,trial_grad,trial_hess)
            if (trial_h >= h .or. abs(step) < 1.0e-10_dp) exit
            step = 0.5_dp*step
            trial = mode+step
         end do
         mode = trial
         h = trial_h
         grad = trial_grad
         hess = trial_hess
         if (abs(step) < 1.0e-9_dp*(1.0_dp+abs(mode))) exit
      end do
      if (hess >= -1.0e-12_dp) then
         info = 1
         return
      end if
      scale = 1.0_dp/sqrt(-hess)
      allocate(log_terms(size(aghq_rule%nodes)))
      do j = 1, size(log_terms)
         call cluster_derivatives(beta,sigma,level, &
            mode+scale*aghq_rule%nodes(j),trial_h,trial_grad,trial_hess)
         log_terms(j) = log(max(tiny(1.0_dp),aghq_rule%weights(j)))+trial_h + &
            0.5_dp*aghq_rule%nodes(j)**2
      end do
      maximum = maxval(log_terms)
      log_integral = log(sqrt(2.0_dp*pi)*scale)+maximum + &
         log(sum(exp(log_terms-maximum)))
   end subroutine cluster_mode_integral

   subroutine cluster_derivatives(beta,sigma,level,b,h,grad,hess)
      real(dp), intent(in) :: beta(:), sigma, b
      integer, intent(in) :: level
      real(dp), intent(out) :: h, grad, hess
      real(dp) :: eta, ll, score, curvature
      integer :: i
      h = -0.5_dp*(b/sigma)**2-log(sigma)-0.5_dp*log(2.0_dp*pi)
      grad = -b/(sigma*sigma)
      hess = -1.0_dp/(sigma*sigma)
      do i = 1, size(aghq_y)
         if (aghq_group(i) /= level) cycle
         eta = aghq_offset(i)+dot_product(aghq_x(i,:),beta)+aghq_z(i)*b
         call observation_terms(aghq_family,aghq_y(i),eta,aghq_dispersion, &
            ll,score,curvature)
         h = h+aghq_weights(i)*ll
         grad = grad+aghq_weights(i)*score*aghq_z(i)
         hess = hess+aghq_weights(i)*curvature*aghq_z(i)**2
      end do
   end subroutine cluster_derivatives

   subroutine observation_terms(family,y,eta,dispersion,ll,score,hess)
      integer, intent(in) :: family
      real(dp), intent(in) :: y, eta, dispersion
      real(dp), intent(out) :: ll, score, hess
      real(dp) :: mu
      ll = -huge(1.0_dp)/100.0_dp
      score = 0.0_dp
      hess = -1.0_dp
      select case (family)
      case (family_binomial)
         mu = min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,logistic(eta)))
         ll = y*log(mu)+(1.0_dp-y)*log(1.0_dp-mu)
         score = y-mu
         hess = -mu*(1.0_dp-mu)
      case (family_poisson)
         mu = exp(min(30.0_dp,max(-30.0_dp,eta)))
         ll = y*log(mu)-mu-log_gamma(y+1.0_dp)
         score = y-mu
         hess = -mu
      case (family_negative_binomial)
         mu = exp(min(30.0_dp,max(-30.0_dp,eta)))
         ll = log_gamma(y+dispersion)-log_gamma(dispersion)-log_gamma(y+1.0_dp) + &
            dispersion*log(dispersion/(dispersion+mu))+y*log(mu/(dispersion+mu))
         score = dispersion*(y-mu)/(dispersion+mu)
         hess = -dispersion*mu*(dispersion+y)/(dispersion+mu)**2
      end select
   end subroutine observation_terms

   subroutine numerical_hessian(par,hessian)
      real(dp), intent(in) :: par(:)
      real(dp), allocatable, intent(out) :: hessian(:,:)
      real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:)
      real(dp) :: f0, fp, fm, fpp, fpm, fmp, fmm, hi, hj
      integer :: n, i, j
      n = size(par)
      allocate(hessian(n,n),xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n))
      hessian = 0.0_dp
      f0 = aghq_objective(par)
      do i = 1, n
         hi = 1.0e-4_dp*(1.0_dp+abs(par(i)))
         xp = par
         xm = par
         xp(i) = xp(i)+hi
         xm(i) = xm(i)-hi
         fp = aghq_objective(xp)
         fm = aghq_objective(xm)
         hessian(i,i) = (fp-2.0_dp*f0+fm)/(hi*hi)
         do j = i+1, n
            hj = 1.0e-4_dp*(1.0_dp+abs(par(j)))
            xpp = par
            xpm = par
            xmp = par
            xmm = par
            xpp(i)=xpp(i)+hi; xpp(j)=xpp(j)+hj
            xpm(i)=xpm(i)+hi; xpm(j)=xpm(j)-hj
            xmp(i)=xmp(i)-hi; xmp(j)=xmp(j)+hj
            xmm(i)=xmm(i)-hi; xmm(j)=xmm(j)-hj
            fpp=aghq_objective(xpp)
            fpm=aghq_objective(xpm)
            fmp=aghq_objective(xmp)
            fmm=aghq_objective(xmm)
            hessian(i,j)=(fpp-fpm-fmp+fmm)/(4.0_dp*hi*hj)
            hessian(j,i)=hessian(i,j)
         end do
      end do
      do i = 1, n
         hessian(i,i) = hessian(i,i)+1.0e-8_dp
      end do
   end subroutine numerical_hessian

   pure real(dp) function logistic(x) result(value)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         value = 1.0_dp/(1.0_dp+exp(-min(700.0_dp,x)))
      else
         value = exp(max(-700.0_dp,x))/(1.0_dp+exp(max(-700.0_dp,x)))
      end if
   end function logistic

   subroutine initialize_result(result,family,dispersion,order)
      type(glmm_result_t), intent(out) :: result
      integer, intent(in) :: family, order
      real(dp), intent(in) :: dispersion
      result%family = family
      result%dispersion = dispersion
      result%quadrature_order = order
      result%message = 'not fitted'
   end subroutine initialize_result

   subroutine fail_result(result,status,message)
      type(glmm_result_t), intent(inout) :: result
      integer, intent(in) :: status
      character(len=*), intent(in) :: message
      result%status = status
      result%converged = .false.
      result%message = message
   end subroutine fail_result

   subroutine clear_active()
      if (allocated(aghq_y)) deallocate(aghq_y)
      if (allocated(aghq_x)) deallocate(aghq_x)
      if (allocated(aghq_z)) deallocate(aghq_z)
      if (allocated(aghq_weights)) deallocate(aghq_weights)
      if (allocated(aghq_offset)) deallocate(aghq_offset)
      if (allocated(aghq_group)) deallocate(aghq_group)
      if (allocated(aghq_rule%nodes)) deallocate(aghq_rule%nodes)
      if (allocated(aghq_rule%weights)) deallocate(aghq_rule%weights)
      if (allocated(aghq_rule%log_density)) deallocate(aghq_rule%log_density)
      aghq_levels = 0
      aghq_dispersion = 1.0_dp
   end subroutine clear_active

end module lme4_aghq
