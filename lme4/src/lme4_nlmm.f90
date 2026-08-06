module lme4_nlmm
   use lme4_kinds, only : dp, pi
   use lme4_types, only : random_term_t, nlmm_control_t, nlmm_result_t, &
      covariance_unstructured, covariance_diagonal, covariance_compound_symmetry, &
      covariance_ar1
   use lme4_covariance, only : term_covariance_from_eta, eta_to_theta
   use lme4_linalg, only : cholesky_lower, chol_solve, invert_spd, logdet_from_chol
   use minqa_module, only : bobyqa, minqa_control_t, minqa_result_t
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   abstract interface
      function nonlinear_mean_function(covariates,beta,random_effect) result(mean)
         import :: dp
         real(dp), intent(in) :: covariates(:), beta(:), random_effect(:)
         real(dp) :: mean
      end function nonlinear_mean_function
   end interface

   public :: nonlinear_mean_function
   public :: fit_nlmm, predict_nlmm, simulate_nlmm

   real(dp), allocatable :: active_y(:), active_covariates(:,:), active_weights(:)
   integer, allocatable :: active_group(:)
   type(random_term_t) :: active_term
   type(nlmm_control_t) :: active_control
   procedure(nonlinear_mean_function), pointer :: active_mean => null()
   integer :: active_levels = 0
   integer :: active_q = 0
   integer :: active_p = 0

contains

   subroutine fit_nlmm(y,covariates,group,n_levels,n_random,mean_function,start_beta, &
      result,weights,control,covariance_structure,start_covariance,sigma_start)
      real(dp), intent(in) :: y(:), covariates(:,:)
      integer, intent(in) :: group(:), n_levels, n_random
      procedure(nonlinear_mean_function) :: mean_function
      real(dp), intent(in) :: start_beta(:)
      type(nlmm_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), start_covariance(:,:), sigma_start
      type(nlmm_control_t), intent(in), optional :: control
      integer, intent(in), optional :: covariance_structure

      type(nlmm_control_t) :: ctrl
      type(minqa_control_t) :: mctrl
      type(minqa_result_t) :: mresult
      type(random_term_t) :: terms(1)
      real(dp), allocatable :: w(:), par(:), lower(:), upper(:), modes(:,:)
      real(dp), allocatable :: covariance(:,:), fitted(:), hessian(:,:), hinv(:,:)
      real(dp) :: criterion, start_sigma
      integer :: n, p, q, nt, npar, info, i, hinfo

      ctrl = nlmm_control_t()
      if (present(control)) ctrl = control
      result%message = 'not fitted'
      n = size(y)
      p = size(start_beta)
      q = n_random
      if (size(covariates,1) /= n .or. size(group) /= n .or. n <= p .or. &
          p < 1 .or. q < 1 .or. n_levels < 1) then
         call fail_result(result,1,'invalid nonlinear mixed-model dimensions')
         return
      end if
      if (any(group < 1) .or. any(group > n_levels)) then
         call fail_result(result,1,'group indices must lie between 1 and n_levels')
         return
      end if
      if (.not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(covariates)) .or. &
          .not. all(ieee_is_finite(start_beta))) then
         call fail_result(result,1,'response, covariates, and starting values must be finite')
         return
      end if
      allocate(w(n))
      w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights <= 0.0_dp)) then
            call fail_result(result,1,'weights must be positive and match y')
            return
         end if
         w = weights
      end if
      allocate(active_term%z(n,q),active_term%group(n))
      active_term%z = 0.0_dp
      active_term%group = group
      active_term%n_levels = n_levels
      active_term%name = 'group'
      active_term%covariance_structure = covariance_unstructured
      if (present(covariance_structure)) active_term%covariance_structure = covariance_structure
      if ((active_term%covariance_structure == covariance_compound_symmetry .or. &
           active_term%covariance_structure == covariance_ar1) .and. q < 2) then
         call fail_result(result,1,'compound-symmetry and AR(1) covariance require n_random >= 2')
         return
      end if
      nt = active_term%n_parameters()
      npar = p+nt+1
      allocate(par(npar),lower(npar),upper(npar))
      par(1:p) = start_beta
      lower(1:p) = ctrl%lower_beta
      upper(1:p) = ctrl%upper_beta
      call initialize_covariance(active_term,par(p+1:p+nt),lower(p+1:p+nt), &
         upper(p+1:p+nt),ctrl,start_covariance)
      start_sigma = sqrt(sum(w*(y-sum(w*y)/sum(w))**2)/sum(w))
      if (present(sigma_start)) start_sigma = sigma_start
      par(npar) = log(max(1.0e-6_dp,start_sigma))
      lower(npar) = ctrl%lower_log_sigma
      upper(npar) = ctrl%upper_log_sigma
      par = min(upper,max(lower,par))

      active_y = y
      active_covariates = covariates
      active_group = group
      active_weights = w
      active_control = ctrl
      active_mean => mean_function
      active_levels = n_levels
      active_q = q
      active_p = p

      mctrl%maxfun = ctrl%maxfun
      mctrl%rhoend = max(1.0e-8_dp,ctrl%tolerance)
      mctrl%rhobeg = min(0.25_dp,0.2_dp*minval(upper-lower))
      mctrl%npt = min((npar+1)*(npar+2)/2,max(npar+2,2*npar+1))
      call bobyqa(nlmm_objective,par,mresult,lower,upper,mctrl)
      criterion = mresult%fval
      if (.not. ieee_is_finite(criterion)) then
         call fail_result(result,2,'nonlinear mixed-model optimization failed')
         call clear_active()
         return
      end if
      allocate(modes(q,n_levels))
      call evaluate_nlmm(par,modes,criterion,info)
      if (info /= 0) then
         call fail_result(result,2,'nonlinear mixed-model mode evaluation failed')
         call clear_active()
         return
      end if
      call term_covariance_from_eta(active_term,par(p+1:p+nt),covariance,info)
      if (info /= 0) then
         call fail_result(result,2,'failed to reconstruct nonlinear random-effect covariance')
         call clear_active()
         return
      end if
      allocate(fitted(n))
      do i = 1, n
         fitted(i) = active_mean(covariates(i,:),par(1:p),modes(:,group(i)))
      end do
      result%beta = par(1:p)
      result%u = modes
      terms(1) = active_term
      call eta_to_theta(terms,par(p+1:p+nt),result%theta)
      result%covariance = covariance
      result%sigma = exp(par(npar))
      result%fitted = fitted
      result%residuals = y-fitted
      result%deviance = criterion
      result%log_likelihood = -0.5_dp*criterion
      result%aic = criterion+2.0_dp*real(npar,dp)
      result%bic = criterion+log(real(n,dp))*real(npar,dp)
      result%evaluations = mresult%evaluations
      result%status = 0
      result%converged = .true.
      result%message = 'converged with Gaussian Laplace nonlinear mixed-model fit'
      call numerical_hessian(par,hessian)
      call invert_spd(hessian,hinv,hinfo)
      allocate(result%vcov_beta(p,p))
      if (hinfo == 0) then
         result%vcov_beta = 2.0_dp*hinv(1:p,1:p)
      else
         result%vcov_beta = 0.0_dp
      end if
      call clear_active()
   end subroutine fit_nlmm

   real(dp) function nlmm_objective(par) result(value)
      real(dp), intent(in) :: par(:)
      real(dp), allocatable :: modes(:,:)
      integer :: info
      allocate(modes(active_q,active_levels))
      call evaluate_nlmm(par,modes,value,info)
      if (info /= 0 .or. .not. ieee_is_finite(value)) value = huge(1.0_dp)/100.0_dp
   end function nlmm_objective

   subroutine evaluate_nlmm(par,modes,criterion,info)
      real(dp), intent(in) :: par(:)
      real(dp), intent(out) :: modes(:,:), criterion
      integer, intent(out) :: info
      real(dp), allocatable :: covariance(:,:), precision(:,:)
      real(dp) :: logdet_cov, sigma, log_integral
      integer :: g, nt

      nt = active_term%n_parameters()
      sigma = exp(par(size(par)))
      call term_covariance_from_eta(active_term,par(active_p+1:active_p+nt),covariance,info)
      if (info /= 0) return
      call invert_spd(covariance,precision,info,logdet_cov)
      if (info /= 0) return
      criterion = 0.0_dp
      log_integral = 0.0_dp
      do g = 1, active_levels
         log_integral = 0.0_dp
         call group_laplace(par(1:active_p),precision,logdet_cov,sigma,g, &
            modes(:,g),log_integral,info)
         if (info /= 0) return
         criterion = criterion-2.0_dp*log_integral
      end do
   end subroutine evaluate_nlmm

   subroutine group_laplace(beta,precision,logdet_cov,sigma,level,mode,log_integral,info)
      real(dp), intent(in) :: beta(:), precision(:,:), logdet_cov, sigma
      integer, intent(in) :: level
      real(dp), intent(out) :: mode(:), log_integral
      integer, intent(out) :: info
      real(dp), allocatable :: gradient(:), hessian(:,:), negative_hessian(:,:), chol(:,:), step(:), trial(:)
      real(dp), allocatable :: trial_gradient(:), trial_hessian(:,:)
      real(dp) :: h, trial_h
      integer :: iter, shrink, q

      log_integral = -huge(1.0_dp)
      q = size(mode)
      mode = 0.0_dp
      call numerical_mode_derivatives(beta,precision,logdet_cov,sigma,level,mode,h,gradient,hessian)
      info = 0
      do iter = 1, active_control%max_mode_iterations
         negative_hessian = -hessian
         call cholesky_lower(negative_hessian,chol,info,jitter=1.0e-8_dp)
         if (info /= 0) return
         call chol_solve(chol,gradient,step)
         if (sqrt(dot_product(step,step)) > max(1.0_dp,2.0_dp*sqrt(real(q,dp)))) then
            step = step*max(1.0_dp,2.0_dp*sqrt(real(q,dp)))/sqrt(dot_product(step,step))
         end if
         trial = mode+step
         do shrink = 1, 20
            call numerical_mode_derivatives(beta,precision,logdet_cov,sigma,level,trial, &
               trial_h,trial_gradient,trial_hessian)
            if (trial_h >= h .or. maxval(abs(step)) < 1.0e-10_dp) exit
            step = 0.5_dp*step
            trial = mode+step
         end do
         mode = trial
         h = trial_h
         gradient = trial_gradient
         hessian = trial_hessian
         if (maxval(abs(step)) < active_control%mode_tolerance*(1.0_dp+maxval(abs(mode)))) exit
      end do
      negative_hessian = -hessian
      call cholesky_lower(negative_hessian,chol,info,jitter=1.0e-8_dp)
      if (info /= 0) return
      log_integral = h+0.5_dp*real(q,dp)*log(2.0_dp*pi)-0.5_dp*logdet_from_chol(chol)
   end subroutine group_laplace

   subroutine numerical_mode_derivatives(beta,precision,logdet_cov,sigma,level,b,h,gradient,hessian)
      real(dp), intent(in) :: beta(:), precision(:,:), logdet_cov, sigma, b(:)
      integer, intent(in) :: level
      real(dp), intent(out) :: h
      real(dp), allocatable, intent(out) :: gradient(:), hessian(:,:)
      real(dp), allocatable :: bp(:), bm(:), bpp(:), bpm(:), bmp(:), bmm(:)
      real(dp) :: fp, fm, fpp, fpm, fmp, fmm, hi, hj
      integer :: q, i, j

      q = size(b)
      allocate(gradient(q),hessian(q,q),bp(q),bm(q),bpp(q),bpm(q),bmp(q),bmm(q))
      call group_log_density(beta,precision,logdet_cov,sigma,level,b,h)
      do i = 1, q
         hi = 1.0e-4_dp*(1.0_dp+abs(b(i)))
         bp=b; bm=b; bp(i)=bp(i)+hi; bm(i)=bm(i)-hi
         call group_log_density(beta,precision,logdet_cov,sigma,level,bp,fp)
         call group_log_density(beta,precision,logdet_cov,sigma,level,bm,fm)
         gradient(i) = (fp-fm)/(2.0_dp*hi)
         hessian(i,i) = (fp-2.0_dp*h+fm)/(hi*hi)
         do j = i+1, q
            hj = 1.0e-4_dp*(1.0_dp+abs(b(j)))
            bpp=b; bpm=b; bmp=b; bmm=b
            bpp(i)=bpp(i)+hi; bpp(j)=bpp(j)+hj
            bpm(i)=bpm(i)+hi; bpm(j)=bpm(j)-hj
            bmp(i)=bmp(i)-hi; bmp(j)=bmp(j)+hj
            bmm(i)=bmm(i)-hi; bmm(j)=bmm(j)-hj
            call group_log_density(beta,precision,logdet_cov,sigma,level,bpp,fpp)
            call group_log_density(beta,precision,logdet_cov,sigma,level,bpm,fpm)
            call group_log_density(beta,precision,logdet_cov,sigma,level,bmp,fmp)
            call group_log_density(beta,precision,logdet_cov,sigma,level,bmm,fmm)
            hessian(i,j) = (fpp-fpm-fmp+fmm)/(4.0_dp*hi*hj)
            hessian(j,i) = hessian(i,j)
         end do
      end do
   end subroutine numerical_mode_derivatives

   subroutine group_log_density(beta,precision,logdet_cov,sigma,level,b,h)
      real(dp), intent(in) :: beta(:), precision(:,:), logdet_cov, sigma, b(:)
      integer, intent(in) :: level
      real(dp), intent(out) :: h
      real(dp) :: mean, residual, variance
      integer :: i, q

      q = size(b)
      h = -0.5_dp*dot_product(b,matmul(precision,b))-0.5_dp*logdet_cov- &
         0.5_dp*real(q,dp)*log(2.0_dp*pi)
      do i = 1, size(active_y)
         if (active_group(i) /= level) cycle
         mean = active_mean(active_covariates(i,:),beta,b)
         if (.not. ieee_is_finite(mean)) then
            h = -huge(1.0_dp)/100.0_dp
            return
         end if
         residual = active_y(i)-mean
         variance = sigma*sigma/active_weights(i)
         h = h-0.5_dp*(log(2.0_dp*pi*variance)+residual*residual/variance)
      end do
   end subroutine group_log_density

   subroutine initialize_covariance(term,eta,lower,upper,ctrl,start_covariance)
      type(random_term_t), intent(in) :: term
      real(dp), intent(out) :: eta(:), lower(:), upper(:)
      type(nlmm_control_t), intent(in) :: ctrl
      real(dp), intent(in), optional :: start_covariance(:,:)
      real(dp), allocatable :: chol(:,:)
      real(dp) :: sigma, rho, lower_rho, probability, shift
      integer :: q, i, j, idx, info

      q = term%n_coefficients()
      idx = 0
      select case (term%covariance_structure)
      case (covariance_unstructured)
         if (present(start_covariance)) then
            call cholesky_lower(start_covariance,chol,info,jitter=1.0e-10_dp)
         else
            info = 1
         end if
         do j=1,q
            do i=j,q
               idx=idx+1
               if (i==j) then
                  eta(idx)=log(0.5_dp)
                  if (info==0) eta(idx)=log(max(1.0e-8_dp,chol(i,j)))
                  lower(idx)=ctrl%lower_log_sd; upper(idx)=ctrl%upper_log_sd
               else
                  eta(idx)=0.0_dp
                  if (info==0) eta(idx)=chol(i,j)
                  lower(idx)=ctrl%lower_offdiag; upper(idx)=ctrl%upper_offdiag
               end if
            end do
         end do
      case (covariance_diagonal)
         do i=1,q
            eta(i)=log(0.5_dp)
            if (present(start_covariance)) eta(i)=0.5_dp*log(max(1.0e-16_dp,start_covariance(i,i)))
            lower(i)=ctrl%lower_log_sd; upper(i)=ctrl%upper_log_sd
         end do
      case (covariance_compound_symmetry)
         sigma=0.5_dp; rho=0.0_dp
         if (present(start_covariance)) then
            sigma=sqrt(max(1.0e-16_dp,sum([(start_covariance(i,i),i=1,q)])/real(q,dp)))
            rho=sum(start_covariance)-sum([(start_covariance(i,i),i=1,q)])
            rho=rho/(real(q*(q-1),dp)*sigma*sigma)
         end if
         eta(1)=log(sigma)
         lower_rho=-1.0_dp/real(q-1,dp)+100.0_dp*epsilon(1.0_dp)
         probability=-lower_rho/(1.0_dp-lower_rho)
         shift=log(probability/(1.0_dp-probability))
         probability=(min(1.0_dp-1.0e-8_dp,max(lower_rho+1.0e-8_dp,rho))-lower_rho)/(1.0_dp-lower_rho)
         eta(2)=log(probability/(1.0_dp-probability))-shift
         lower=[ctrl%lower_log_sd,ctrl%lower_offdiag]
         upper=[ctrl%upper_log_sd,ctrl%upper_offdiag]
      case (covariance_ar1)
         sigma=0.5_dp; rho=0.0_dp
         if (present(start_covariance)) then
            sigma=sqrt(max(1.0e-16_dp,sum([(start_covariance(i,i),i=1,q)])/real(q,dp)))
            rho=sum([(start_covariance(i,i+1),i=1,q-1)])/real(q-1,dp)/(sigma*sigma)
         end if
         eta=[log(sigma),atanh(min(0.999999_dp,max(-0.999999_dp,rho)))]
         lower=[ctrl%lower_log_sd,ctrl%lower_offdiag]
         upper=[ctrl%upper_log_sd,ctrl%upper_offdiag]
      end select
   end subroutine initialize_covariance

   subroutine numerical_hessian(par,hessian)
      real(dp), intent(in) :: par(:)
      real(dp), allocatable, intent(out) :: hessian(:,:)
      real(dp), allocatable :: xp(:),xm(:),xpp(:),xpm(:),xmp(:),xmm(:)
      real(dp) :: f0,fp,fm,fpp,fpm,fmp,fmm,hi,hj
      integer :: n,i,j
      n=size(par)
      allocate(hessian(n,n),xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n))
      hessian=0.0_dp
      f0=nlmm_objective(par)
      do i=1,n
         hi=1.0e-4_dp*(1.0_dp+abs(par(i)))
         xp=par; xm=par; xp(i)=xp(i)+hi; xm(i)=xm(i)-hi
         fp=nlmm_objective(xp); fm=nlmm_objective(xm)
         hessian(i,i)=(fp-2.0_dp*f0+fm)/(hi*hi)
         do j=i+1,n
            hj=1.0e-4_dp*(1.0_dp+abs(par(j)))
            xpp=par; xpm=par; xmp=par; xmm=par
            xpp(i)=xpp(i)+hi; xpp(j)=xpp(j)+hj
            xpm(i)=xpm(i)+hi; xpm(j)=xpm(j)-hj
            xmp(i)=xmp(i)-hi; xmp(j)=xmp(j)+hj
            xmm(i)=xmm(i)-hi; xmm(j)=xmm(j)-hj
            fpp=nlmm_objective(xpp); fpm=nlmm_objective(xpm)
            fmp=nlmm_objective(xmp); fmm=nlmm_objective(xmm)
            hessian(i,j)=(fpp-fpm-fmp+fmm)/(4.0_dp*hi*hj)
            hessian(j,i)=hessian(i,j)
         end do
      end do
      do i=1,n
         hessian(i,i)=hessian(i,i)+1.0e-8_dp
      end do
   end subroutine numerical_hessian

   subroutine predict_nlmm(result,covariates,group,mean_function,prediction,include_random)
      type(nlmm_result_t), intent(in) :: result
      real(dp), intent(in) :: covariates(:,:)
      integer, intent(in) :: group(:)
      procedure(nonlinear_mean_function) :: mean_function
      real(dp), allocatable, intent(out) :: prediction(:)
      logical, intent(in), optional :: include_random
      logical :: use_random
      real(dp), allocatable :: zero(:)
      integer :: i
      use_random=.true.
      if (present(include_random)) use_random=include_random
      allocate(prediction(size(group)),zero(size(result%u,1)),source=0.0_dp)
      do i=1,size(group)
         if (use_random) then
            prediction(i)=mean_function(covariates(i,:),result%beta,result%u(:,group(i)))
         else
            prediction(i)=mean_function(covariates(i,:),result%beta,zero)
         end if
      end do
   end subroutine predict_nlmm

   subroutine simulate_nlmm(covariates,group,beta,covariance,sigma,mean_function,y, &
      random_effects,seed)
      real(dp), intent(in) :: covariates(:,:), beta(:), covariance(:,:), sigma
      integer, intent(in) :: group(:)
      procedure(nonlinear_mean_function) :: mean_function
      real(dp), allocatable, intent(out) :: y(:), random_effects(:,:)
      integer, intent(in), optional :: seed
      real(dp), allocatable :: chol(:,:), z(:)
      integer :: n_levels, q, info, g, i
      if (present(seed)) call set_seed(seed)
      n_levels=maxval(group)
      q=size(covariance,1)
      call cholesky_lower(covariance,chol,info,jitter=1.0e-10_dp)
      allocate(random_effects(q,n_levels),z(q),y(size(group)))
      do g=1,n_levels
         call normal_vector(z)
         random_effects(:,g)=matmul(chol,z)
      end do
      do i=1,size(group)
         call normal_vector(z(1:1))
         y(i)=mean_function(covariates(i,:),beta,random_effects(:,group(i)))+sigma*z(1)
      end do
   end subroutine simulate_nlmm

   subroutine normal_vector(x)
      real(dp), intent(out) :: x(:)
      real(dp) :: u1,u2,r,angle
      integer :: i
      i=1
      do while (i<=size(x))
         call random_number(u1); call random_number(u2)
         u1=max(tiny(1.0_dp),u1)
         r=sqrt(-2.0_dp*log(u1)); angle=2.0_dp*pi*u2
         x(i)=r*cos(angle)
         if (i+1<=size(x)) x(i+1)=r*sin(angle)
         i=i+2
      end do
   end subroutine normal_vector

   subroutine set_seed(seed)
      integer, intent(in) :: seed
      integer, allocatable :: values(:)
      integer :: n,i
      call random_seed(size=n)
      allocate(values(n))
      do i=1,n
         values(i)=mod(abs(seed)+104729*i,huge(1)-1)+1
      end do
      call random_seed(put=values)
   end subroutine set_seed

   subroutine fail_result(result,status,message)
      type(nlmm_result_t), intent(inout) :: result
      integer, intent(in) :: status
      character(len=*), intent(in) :: message
      result%status=status
      result%converged=.false.
      result%message=message
   end subroutine fail_result

   subroutine clear_active()
      if (allocated(active_y)) deallocate(active_y)
      if (allocated(active_covariates)) deallocate(active_covariates)
      if (allocated(active_weights)) deallocate(active_weights)
      if (allocated(active_group)) deallocate(active_group)
      if (allocated(active_term%z)) deallocate(active_term%z)
      if (allocated(active_term%group)) deallocate(active_term%group)
      if (allocated(active_term%name)) deallocate(active_term%name)
      nullify(active_mean)
      active_levels=0
      active_q=0
      active_p=0
   end subroutine clear_active

end module lme4_nlmm
