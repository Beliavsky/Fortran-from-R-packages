module lme4_aghq_nd
   use lme4_kinds, only : dp, pi
   use lme4_types, only : random_term_t, covariance_block_t, glmm_control_t, &
      glmm_result_t, gh_rule_t, family_binomial, family_poisson, &
      family_negative_binomial, covariance_unstructured, covariance_diagonal, &
      covariance_compound_symmetry, covariance_ar1
   use lme4_covariance, only : validate_terms, term_covariance_from_eta, &
      eta_to_theta, cov2sdcor
   use lme4_quadrature, only : gh_rule
   use lme4_linalg, only : cholesky_lower, chol_solve, invert_spd, logdet_from_chol
   use lme4_glmm, only : fit_glmm
   use minqa_module, only : bobyqa, minqa_control_t, minqa_result_t
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private
   public :: fit_glmm_aghq_multidimensional

   real(dp), allocatable :: active_y(:), active_x(:,:), active_z(:,:)
   real(dp), allocatable :: active_weights(:), active_offset(:)
   integer, allocatable :: active_group(:)
   type(random_term_t) :: active_term
   type(gh_rule_t) :: active_rule
   integer :: active_family = 0
   integer :: active_levels = 0
   integer :: active_q = 0
   real(dp) :: active_dispersion = 1.0_dp

contains

   subroutine fit_glmm_aghq_multidimensional(y,x,term,family,result,order,weights, &
      offset,dispersion,control,max_nodes)
      real(dp), intent(in) :: y(:), x(:,:)
      type(random_term_t), intent(inout) :: term
      integer, intent(in) :: family, order
      type(glmm_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), offset(:), dispersion
      type(glmm_control_t), intent(in), optional :: control
      integer, intent(in), optional :: max_nodes

      type(glmm_control_t) :: ctrl
      type(glmm_result_t) :: initial
      type(minqa_control_t) :: mctrl
      type(minqa_result_t) :: mresult
      type(random_term_t) :: terms(1)
      real(dp), allocatable :: w(:), off(:), par(:), lower(:), upper(:)
      real(dp), allocatable :: modes(:,:), eta(:), mu(:), covariance(:,:), sdcor(:,:)
      real(dp), allocatable :: hessian(:,:), hinv(:,:)
      real(dp) :: disp, criterion
      integer :: n, p, q, nt, npar, info, i, g, hinfo, node_limit
      integer(kind=8) :: node_count
      logical :: ok
      character(len=:), allocatable :: message

      ctrl = glmm_control_t()
      if (present(control)) ctrl = control
      disp = 1.0_dp
      if (present(dispersion)) disp = dispersion
      result%family = family
      result%dispersion = disp
      result%quadrature_order = order
      result%message = 'not fitted'
      n = size(y)
      p = size(x,2)
      q = term%n_coefficients()
      node_limit = 20000
      if (present(max_nodes)) node_limit = max_nodes
      if (order < 1 .or. q < 1) then
         call fail_result(result,1,'quadrature order and random-effect dimension must be positive')
         return
      end if
      node_count = int(order,8)**q
      if (node_count > int(node_limit,8)) then
         call fail_result(result,1,'tensor quadrature exceeds max_nodes; reduce order or random-effect dimension')
         return
      end if
      if (size(x,1) /= n .or. p < 1 .or. n <= p) then
         call fail_result(result,1,'invalid fixed-effect design matrix')
         return
      end if
      if (family /= family_binomial .and. family /= family_poisson .and. &
          family /= family_negative_binomial) then
         call fail_result(result,1,'multidimensional AGHQ supports binomial, poisson, and negative binomial')
         return
      end if
      if (family == family_binomial .and. (any(y < 0.0_dp) .or. any(y > 1.0_dp))) then
         call fail_result(result,1,'binomial responses must lie in [0,1]')
         return
      end if
      if ((family == family_poisson .or. family == family_negative_binomial) .and. any(y < 0.0_dp)) then
         call fail_result(result,1,'count responses must be nonnegative')
         return
      end if
      if (family == family_negative_binomial .and. disp <= 0.0_dp) then
         call fail_result(result,1,'negative-binomial size must be positive')
         return
      end if
      terms(1) = term
      call validate_terms(terms,n,ok,message)
      if (.not. ok) then
         call fail_result(result,1,message)
         return
      end if
      term = terms(1)
      allocate(w(n),off(n))
      w = 1.0_dp
      off = 0.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights <= 0.0_dp)) then
            call fail_result(result,1,'weights must be positive and match y')
            return
         end if
         w = weights
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            call fail_result(result,1,'offset must match y')
            return
         end if
         off = offset
      end if
      call gh_rule(order,active_rule,info)
      if (info /= 0) then
         call fail_result(result,1,'failed to construct Gauss-Hermite rule')
         return
      end if

      call fit_glmm(y,x,terms,family,initial,weights=w,offset=off,control=ctrl,dispersion=disp)
      nt = term%n_parameters()
      npar = p+nt
      allocate(par(npar),lower(npar),upper(npar))
      par = 0.0_dp
      if (initial%converged) par(1:p) = initial%beta
      if (.not. initial%converged) then
         if (family == family_binomial) then
            par(1) = log(max(1.0e-6_dp,min(1.0_dp-1.0e-6_dp,sum(w*y)/sum(w)))/ &
               max(1.0e-6_dp,1.0_dp-sum(w*y)/sum(w)))
         else
            par(1) = log(max(1.0e-6_dp,sum(w*y)/sum(w)))
         end if
      end if
      lower(1:p) = ctrl%lower_beta
      upper(1:p) = ctrl%upper_beta
      call initialize_covariance_parameters(term,par(p+1:npar),lower(p+1:npar), &
         upper(p+1:npar),ctrl)

      active_y = y
      active_x = x
      active_z = term%z
      active_group = term%group
      active_weights = w
      active_offset = off
      active_term = term
      active_family = family
      active_dispersion = disp
      active_levels = term%n_levels
      active_q = q

      mctrl%maxfun = ctrl%maxfun
      mctrl%rhoend = max(1.0e-8_dp,ctrl%tolerance)
      mctrl%rhobeg = min(0.25_dp,0.2_dp*minval(upper-lower))
      mctrl%npt = min((npar+1)*(npar+2)/2,max(npar+2,2*npar+1))
      call bobyqa(nd_objective,par,mresult,lower,upper,mctrl)
      criterion = mresult%fval
      if (.not. ieee_is_finite(criterion)) then
         call fail_result(result,2,'multidimensional AGHQ optimization failed')
         call clear_active()
         return
      end if
      allocate(modes(q,active_levels))
      call evaluate_all_clusters(par(1:p),par(p+1:npar),modes,criterion,info)
      if (info /= 0) then
         call fail_result(result,2,'failed to evaluate multidimensional adaptive quadrature')
         call clear_active()
         return
      end if
      call term_covariance_from_eta(term,par(p+1:npar),covariance,info)
      if (info /= 0) then
         call fail_result(result,2,'failed to reconstruct random-effect covariance')
         call clear_active()
         return
      end if
      call cov2sdcor(covariance,sdcor)

      result%beta = par(1:p)
      allocate(result%u(q*active_levels))
      do g = 1, active_levels
         result%u((g-1)*q+1:g*q) = modes(:,g)
      end do
      terms(1) = term
      call eta_to_theta(terms,par(p+1:npar),result%theta)
      allocate(result%varcorr(1))
      result%varcorr(1)%name = term%name
      result%varcorr(1)%n_levels = term%n_levels
      result%varcorr(1)%covariance = covariance
      result%varcorr(1)%sdcor = sdcor
      allocate(result%term_offsets(2))
      result%term_offsets = [1,q*active_levels+1]
      eta = off+matmul(x,result%beta)
      do i = 1, n
         eta(i) = eta(i)+dot_product(term%z(i,:),modes(:,term%group(i)))
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
      result%evaluations = mresult%evaluations
      result%status = 0
      result%converged = .true.
      result%message = 'converged with multidimensional adaptive Gauss-Hermite quadrature'

      call numerical_hessian(par,hessian)
      call invert_spd(hessian,hinv,hinfo)
      allocate(result%vcov_beta(p,p))
      if (hinfo == 0) then
         result%vcov_beta = 2.0_dp*hinv(1:p,1:p)
      else
         result%vcov_beta = 0.0_dp
      end if
      call clear_active()
   end subroutine fit_glmm_aghq_multidimensional

   real(dp) function nd_objective(par) result(value)
      real(dp), intent(in) :: par(:)
      real(dp), allocatable :: modes(:,:)
      integer :: p, info
      p = size(active_x,2)
      allocate(modes(active_q,active_levels))
      call evaluate_all_clusters(par(1:p),par(p+1:),modes,value,info)
      if (info /= 0 .or. .not. ieee_is_finite(value)) value = huge(1.0_dp)/100.0_dp
   end function nd_objective

   subroutine evaluate_all_clusters(beta,eta_cov,modes,criterion,info)
      real(dp), intent(in) :: beta(:), eta_cov(:)
      real(dp), intent(out) :: modes(:,:), criterion
      integer, intent(out) :: info
      real(dp), allocatable :: covariance(:,:), precision(:,:)
      real(dp) :: logdet_cov, log_integral
      integer :: g

      call term_covariance_from_eta(active_term,eta_cov,covariance,info)
      if (info /= 0) return
      call invert_spd(covariance,precision,info,logdet_cov)
      if (info /= 0) return
      criterion = 0.0_dp
      log_integral = 0.0_dp
      do g = 1, active_levels
         log_integral = 0.0_dp
         call cluster_mode_integral(beta,precision,logdet_cov,g,modes(:,g),log_integral,info)
         if (info /= 0) return
         criterion = criterion-2.0_dp*log_integral
      end do
   end subroutine evaluate_all_clusters

   subroutine cluster_mode_integral(beta,precision,logdet_cov,level,mode,log_integral,info)
      real(dp), intent(in) :: beta(:), precision(:,:), logdet_cov
      integer, intent(in) :: level
      real(dp), intent(out) :: mode(:), log_integral
      integer, intent(out) :: info
      real(dp), allocatable :: gradient(:), hessian(:,:), negative_hessian(:,:), chol(:,:)
      real(dp), allocatable :: step(:), trial(:), trial_gradient(:), trial_hessian(:,:)
      real(dp), allocatable :: inverse_transpose(:,:), identity(:,:), node(:), b(:), log_terms(:)
      real(dp) :: h, trial_h, maximum, log_jacobian
      integer :: iter, shrink, q, n_nodes, index, d, tmp, node_index

      log_integral = -huge(1.0_dp)
      q = size(mode)
      mode = 0.0_dp
      call cluster_derivatives(beta,precision,logdet_cov,level,mode,h,gradient,hessian)
      info = 0
      do iter = 1, 60
         negative_hessian = -hessian
         call cholesky_lower(negative_hessian,chol,info,jitter=1.0e-10_dp)
         if (info /= 0) return
         call chol_solve(chol,gradient,step)
         if (sqrt(dot_product(step,step)) > max(1.0_dp,2.0_dp*sqrt(real(q,dp)))) then
            step = step*max(1.0_dp,2.0_dp*sqrt(real(q,dp)))/sqrt(dot_product(step,step))
         end if
         trial = mode+step
         do shrink = 1, 20
            call cluster_derivatives(beta,precision,logdet_cov,level,trial, &
               trial_h,trial_gradient,trial_hessian)
            if (trial_h >= h .or. maxval(abs(step)) < 1.0e-10_dp) exit
            step = 0.5_dp*step
            trial = mode+step
         end do
         mode = trial
         h = trial_h
         gradient = trial_gradient
         hessian = trial_hessian
         if (maxval(abs(step)) < 1.0e-9_dp*(1.0_dp+maxval(abs(mode)))) exit
      end do
      negative_hessian = -hessian
      call cholesky_lower(negative_hessian,chol,info,jitter=1.0e-10_dp)
      if (info /= 0) return
      allocate(identity(q,q),source=0.0_dp)
      do d = 1, q
         identity(d,d) = 1.0_dp
      end do
      call solve_upper_from_lower(chol,identity,inverse_transpose)
      log_jacobian = -sum(log([(chol(d,d),d=1,q)]))
      n_nodes = size(active_rule%nodes)**q
      allocate(log_terms(n_nodes),node(q),b(q))
      do index = 0, n_nodes-1
         tmp = index
         log_terms(index+1) = 0.0_dp
         do d = 1, q
            node_index = mod(tmp,size(active_rule%nodes))+1
            tmp = tmp/size(active_rule%nodes)
            node(d) = active_rule%nodes(node_index)
            log_terms(index+1) = log_terms(index+1)+ &
               log(max(tiny(1.0_dp),active_rule%weights(node_index)))
         end do
         b = mode+matmul(inverse_transpose,node)
         call cluster_log_density(beta,precision,logdet_cov,level,b,trial_h)
         log_terms(index+1) = log_terms(index+1)+trial_h+0.5_dp*dot_product(node,node)
      end do
      maximum = maxval(log_terms)
      log_integral = 0.5_dp*real(q,dp)*log(2.0_dp*pi)+log_jacobian+maximum+ &
         log(sum(exp(log_terms-maximum)))
   end subroutine cluster_mode_integral

   subroutine cluster_derivatives(beta,precision,logdet_cov,level,b,h,gradient,hessian)
      real(dp), intent(in) :: beta(:), precision(:,:), logdet_cov, b(:)
      integer, intent(in) :: level
      real(dp), intent(out) :: h
      real(dp), allocatable, intent(out) :: gradient(:), hessian(:,:)
      real(dp) :: eta, ll, score, curvature
      integer :: i, q

      q = size(b)
      allocate(gradient(q),hessian(q,q))
      h = -0.5_dp*dot_product(b,matmul(precision,b))-0.5_dp*logdet_cov- &
         0.5_dp*real(q,dp)*log(2.0_dp*pi)
      gradient = -matmul(precision,b)
      hessian = -precision
      do i = 1, size(active_y)
         if (active_group(i) /= level) cycle
         eta = active_offset(i)+dot_product(active_x(i,:),beta)+dot_product(active_z(i,:),b)
         call observation_terms(active_family,active_y(i),eta,active_dispersion,ll,score,curvature)
         h = h+active_weights(i)*ll
         gradient = gradient+active_weights(i)*score*active_z(i,:)
         hessian = hessian+active_weights(i)*curvature*outer(active_z(i,:),active_z(i,:))
      end do
   end subroutine cluster_derivatives

   subroutine cluster_log_density(beta,precision,logdet_cov,level,b,h)
      real(dp), intent(in) :: beta(:), precision(:,:), logdet_cov, b(:)
      integer, intent(in) :: level
      real(dp), intent(out) :: h
      real(dp) :: eta, ll, score, curvature
      integer :: i, q
      q = size(b)
      h = -0.5_dp*dot_product(b,matmul(precision,b))-0.5_dp*logdet_cov- &
         0.5_dp*real(q,dp)*log(2.0_dp*pi)
      do i = 1, size(active_y)
         if (active_group(i) /= level) cycle
         eta = active_offset(i)+dot_product(active_x(i,:),beta)+dot_product(active_z(i,:),b)
         call observation_terms(active_family,active_y(i),eta,active_dispersion,ll,score,curvature)
         h = h+active_weights(i)*ll
      end do
   end subroutine cluster_log_density

   subroutine observation_terms(family,y,eta,dispersion,ll,score,hess)
      integer, intent(in) :: family
      real(dp), intent(in) :: y, eta, dispersion
      real(dp), intent(out) :: ll, score, hess
      real(dp) :: mu
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
         ll = log_gamma(y+dispersion)-log_gamma(dispersion)-log_gamma(y+1.0_dp)+ &
            dispersion*log(dispersion/(dispersion+mu))+y*log(mu/(dispersion+mu))
         score = dispersion*(y-mu)/(dispersion+mu)
         hess = -dispersion*mu*(dispersion+y)/(dispersion+mu)**2
      end select
   end subroutine observation_terms

   subroutine solve_upper_from_lower(lower,b,x)
      real(dp), intent(in) :: lower(:,:), b(:,:)
      real(dp), allocatable, intent(out) :: x(:,:)
      integer :: n, i, j, k
      n = size(lower,1)
      allocate(x(n,size(b,2)))
      do j = 1, size(b,2)
         do i = n, 1, -1
            x(i,j) = b(i,j)
            do k = i+1, n
               x(i,j) = x(i,j)-lower(k,i)*x(k,j)
            end do
            x(i,j) = x(i,j)/lower(i,i)
         end do
      end do
   end subroutine solve_upper_from_lower

   pure function outer(a,b) result(matrix)
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: matrix(size(a),size(b))
      integer :: i, j
      do j = 1, size(b)
         do i = 1, size(a)
            matrix(i,j) = a(i)*b(j)
         end do
      end do
   end function outer

   subroutine initialize_covariance_parameters(term,eta,lower,upper,ctrl)
      type(random_term_t), intent(in) :: term
      real(dp), intent(out) :: eta(:), lower(:), upper(:)
      type(glmm_control_t), intent(in) :: ctrl
      integer :: q, i, j, idx
      q = term%n_coefficients()
      idx = 0
      select case (term%covariance_structure)
      case (covariance_unstructured)
         do j = 1, q
            do i = j, q
               idx = idx+1
               if (i == j) then
                  eta(idx)=log(0.5_dp); lower(idx)=ctrl%lower_log_sd; upper(idx)=ctrl%upper_log_sd
               else
                  eta(idx)=0.0_dp; lower(idx)=ctrl%lower_offdiag; upper(idx)=ctrl%upper_offdiag
               end if
            end do
         end do
      case (covariance_diagonal)
         do i = 1, q
            idx=idx+1; eta(idx)=log(0.5_dp); lower(idx)=ctrl%lower_log_sd; upper(idx)=ctrl%upper_log_sd
         end do
      case (covariance_compound_symmetry,covariance_ar1)
         eta(1)=log(0.5_dp); lower(1)=ctrl%lower_log_sd; upper(1)=ctrl%upper_log_sd
         eta(2)=0.0_dp; lower(2)=ctrl%lower_offdiag; upper(2)=ctrl%upper_offdiag
      end select
   end subroutine initialize_covariance_parameters

   subroutine numerical_hessian(par,hessian)
      real(dp), intent(in) :: par(:)
      real(dp), allocatable, intent(out) :: hessian(:,:)
      real(dp), allocatable :: xp(:),xm(:),xpp(:),xpm(:),xmp(:),xmm(:)
      real(dp) :: f0,fp,fm,fpp,fpm,fmp,fmm,hi,hj
      integer :: n,i,j
      n=size(par)
      allocate(hessian(n,n),xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n))
      hessian=0.0_dp
      f0=nd_objective(par)
      do i=1,n
         hi=1.0e-4_dp*(1.0_dp+abs(par(i)))
         xp=par; xm=par; xp(i)=xp(i)+hi; xm(i)=xm(i)-hi
         fp=nd_objective(xp); fm=nd_objective(xm)
         hessian(i,i)=(fp-2.0_dp*f0+fm)/(hi*hi)
         do j=i+1,n
            hj=1.0e-4_dp*(1.0_dp+abs(par(j)))
            xpp=par; xpm=par; xmp=par; xmm=par
            xpp(i)=xpp(i)+hi; xpp(j)=xpp(j)+hj
            xpm(i)=xpm(i)+hi; xpm(j)=xpm(j)-hj
            xmp(i)=xmp(i)-hi; xmp(j)=xmp(j)+hj
            xmm(i)=xmm(i)-hi; xmm(j)=xmm(j)-hj
            fpp=nd_objective(xpp); fpm=nd_objective(xpm)
            fmp=nd_objective(xmp); fmm=nd_objective(xmm)
            hessian(i,j)=(fpp-fpm-fmp+fmm)/(4.0_dp*hi*hj)
            hessian(j,i)=hessian(i,j)
         end do
      end do
      do i=1,n
         hessian(i,i)=hessian(i,i)+1.0e-8_dp
      end do
   end subroutine numerical_hessian

   pure real(dp) function logistic(x) result(value)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         value=1.0_dp/(1.0_dp+exp(-min(700.0_dp,x)))
      else
         value=exp(max(-700.0_dp,x))/(1.0_dp+exp(max(-700.0_dp,x)))
      end if
   end function logistic

   subroutine fail_result(result,status,message)
      type(glmm_result_t), intent(inout) :: result
      integer, intent(in) :: status
      character(len=*), intent(in) :: message
      result%status=status
      result%converged=.false.
      result%message=message
   end subroutine fail_result

   subroutine clear_active()
      if (allocated(active_y)) deallocate(active_y)
      if (allocated(active_x)) deallocate(active_x)
      if (allocated(active_z)) deallocate(active_z)
      if (allocated(active_weights)) deallocate(active_weights)
      if (allocated(active_offset)) deallocate(active_offset)
      if (allocated(active_group)) deallocate(active_group)
      if (allocated(active_rule%nodes)) deallocate(active_rule%nodes)
      if (allocated(active_rule%weights)) deallocate(active_rule%weights)
      if (allocated(active_rule%log_density)) deallocate(active_rule%log_density)
      active_levels=0
      active_q=0
      active_dispersion=1.0_dp
   end subroutine clear_active

end module lme4_aghq_nd
