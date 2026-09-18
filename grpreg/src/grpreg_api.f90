module grpreg_api
   use grpreg_kinds, only : dp
   use grpreg_types, only : grpreg_fit_type, grpreg_selection_type
   use grpreg_preprocess, only : grpreg_preprocess_type, prepare_design, restore_coefficients, normalize_groups
   use grpreg_fit_core, only : fit_path_gaussian, fit_path_glm, fit_path_cox, make_lambda_sequence
   use grpreg_math, only : logistic
   implicit none
   private
   public :: grpreg_fit, gbridge_fit, grpsurv_fit
   public :: coef_grpreg, predict_grpreg, predict_grpsurv_link
   public :: breslow_baseline, predict_grpsurv_survival, predict_grpsurv_hazard, predict_grpsurv_median
   public :: loglik_grpreg, residuals_grpreg, select_grpreg
   public :: count_nonzero, count_nonzero_groups, group_norms

contains

   subroutine grpreg_fit(x, y, group, fit, family, penalty, lambda, nlambda, lambda_min, alpha, gamma, tau, eps, &
      max_iter, group_multiplier, log_lambda, bridge_delta)
      real(dp), intent(in) :: x(:,:) !! Raw design matrix, observations by predictors, without an intercept.
      real(dp), intent(in) :: y(:) !! Response vector for Gaussian, binomial, or Poisson regression.
      integer, intent(in) :: group(:) !! Predictor group labels; zero identifies unpenalized covariates.
      type(grpreg_fit_type), intent(out) :: fit !! Fitted regularization path on the original predictor scale.
      character(len=*), optional, intent(in) :: family !! Outcome family: gaussian, binomial, or poisson; default gaussian.
      character(len=*), optional, intent(in) :: penalty !! Penalty name; default grLasso.
      real(dp), optional, intent(in) :: lambda(:) !! Optional user-supplied regularization sequence.
      integer, optional, intent(in) :: nlambda !! Number of automatically generated lambda values; default 100.
      real(dp), optional, intent(in) :: lambda_min !! Minimum lambda fraction; default depends on n versus p.
      real(dp), optional, intent(in) :: alpha !! Group-penalty fraction in (0,1]; default one.
      real(dp), optional, intent(in) :: gamma !! MCP/SCAD concavity parameter; defaults to three or four for SCAD.
      real(dp), optional, intent(in) :: tau !! Group exponential lasso decay parameter; default one third.
      real(dp), optional, intent(in) :: eps !! Coordinate convergence tolerance; default 1e-4.
      integer, optional, intent(in) :: max_iter !! Maximum coordinate sweeps over the path; default 10000.
      real(dp), optional, intent(in) :: group_multiplier(:) !! Optional nonnegative multiplier per penalized group.
      logical, optional, intent(in) :: log_lambda !! Whether automatic lambda values are log-spaced; default true.
      real(dp), optional, intent(in) :: bridge_delta !! Positive deletion threshold used only by group bridge; default 1e-7.
      character(len=16) :: fam, pen
      type(grpreg_preprocess_type) :: prep
      real(dp), allocatable :: xw(:,:), yw(:), lam(:), bw(:,:), iw(:), dev(:), dfr(:), etaw(:,:), score(:)
      real(dp), allocatable :: beta(:,:), intercept(:)
      integer, allocatable :: it(:), gnorm(:)
      real(dp) :: a, gam, tv, tol, lmin, ymean, bdelta
      integer :: nl, mit, ng
      logical :: bilevel, logl
      if (size(x,1) /= size(y)) error stop 'grpreg: X and y have incompatible dimensions'
      if (size(x,2) /= size(group)) error stop 'grpreg: group length does not match X columns'
      fam = 'gaussian'
      if (present(family)) fam = trim(family)
      pen = 'grLasso'
      if (present(penalty)) pen = trim(penalty)
      if (.not. valid_family(fam)) error stop 'grpreg: unsupported family'
      if (.not. valid_penalty(pen, .true.)) error stop 'grpreg: unsupported penalty'
      a = 1.0_dp
      if (present(alpha)) a = alpha
      if (a <= 0.0_dp .or. a > 1.0_dp) error stop 'grpreg: alpha must lie in (0,1]'
      gam = merge(4.0_dp,3.0_dp,trim(pen) == 'grSCAD')
      if (present(gamma)) gam = gamma
      tv = 1.0_dp/3.0_dp
      if (present(tau)) tv = tau
      tol = 1.0e-4_dp
      if (present(eps)) tol = eps
      mit = 10000
      if (present(max_iter)) mit = max_iter
      nl = 100
      if (present(nlambda)) nl = nlambda
      lmin = merge(0.0001_dp,0.05_dp,size(x,1) > size(x,2))
      if (present(lambda_min)) lmin = lambda_min
      bdelta = 1.0e-7_dp
      if (present(bridge_delta)) bdelta = bridge_delta
      if (bdelta <= 0.0_dp) error stop 'grpreg: bridge_delta must be positive'
      logl = .true.
      if (present(log_lambda)) logl = log_lambda
      bilevel = pen(1:2) /= 'gr'
      call prepare_design(x, group, bilevel, xw, prep, group_multiplier)
      allocate(yw(size(y)))
      ymean = 0.0_dp
      if (trim(fam) == 'gaussian') then
         ymean = sum(y)/real(size(y),dp)
         yw = y-ymean
      else
         yw = y
      end if
      if (present(lambda)) then
         allocate(lam(size(lambda)))
         lam = lambda
      else
         call null_score_glm(xw,yw,prep%group_reduced,fam,score)
         call make_lambda_sequence(xw,score,prep%group_reduced,prep%group_multiplier,a,bilevel,nl,lmin,logl,lam)
      end if
      if (trim(fam) == 'gaussian') then
         call fit_path_gaussian(xw,yw,prep%group_reduced,prep%group_multiplier,pen,lam,a,gam,tv,bdelta,tol,mit, &
            bw,iw,dev,dfr,it,etaw)
      else
         call fit_path_glm(xw,yw,prep%group_reduced,prep%group_multiplier,fam,pen,lam,a,gam,tv,bdelta,tol,mit, &
            bw,iw,dev,dfr,it,etaw)
      end if
      call restore_coefficients(bw,iw,prep,beta,intercept)
      if (trim(fam) == 'gaussian') then
         intercept = intercept+ymean
         dfr = dfr+1.0_dp
      end if
      call normalize_groups(group,gnorm,ng)
      fit%family = fam
      fit%penalty = pen
      fit%n = size(x,1)
      fit%p = size(x,2)
      fit%ngroups = ng
      fit%nlambda = size(lam)
      fit%alpha = a
      fit%gamma = gam
      fit%tau = tv
      fit%lambda = lam
      fit%intercept = intercept
      fit%beta = beta
      fit%deviance = dev
      fit%df = dfr
      fit%iter = it
      allocate(fit%eta(size(y),size(lam)))
      fit%eta = spread(intercept,1,size(y)) + matmul(x,beta)
      fit%group = gnorm
      fit%group_multiplier = prep%group_multiplier
      fit%x_center = prep%center
      fit%x_scale = prep%scale
      fit%y = y
   end subroutine grpreg_fit

   subroutine gbridge_fit(x, y, group, fit, family, lambda, nlambda, lambda_min, lambda_max, alpha, gamma, delta, &
      eps, max_iter, group_multiplier)
      real(dp), intent(in) :: x(:,:) !! Raw design matrix for group bridge regression.
      real(dp), intent(in) :: y(:) !! Gaussian, binomial, or Poisson response vector.
      integer, intent(in) :: group(:) !! Predictor group labels; zero identifies unpenalized covariates.
      type(grpreg_fit_type), intent(out) :: fit !! Fitted group-bridge regularization path.
      character(len=*), optional, intent(in) :: family !! Outcome family; default gaussian.
      real(dp), optional, intent(in) :: lambda(:) !! Optional user-specified group-bridge lambda sequence.
      integer, optional, intent(in) :: nlambda !! Number of generated lambda values; default 100.
      real(dp), optional, intent(in) :: lambda_min !! Minimum generated lambda fraction; default follows upstream n versus p rule.
      real(dp), optional, intent(in) :: lambda_max !! Optional absolute maximum generated lambda.
      real(dp), optional, intent(in) :: alpha !! Sparsity-versus-ridge mixing proportion; default one.
      real(dp), optional, intent(in) :: gamma !! Group-bridge exponent, normally between zero and one; default 0.5.
      real(dp), optional, intent(in) :: delta !! Positive group deletion tolerance; default 1e-7.
      real(dp), optional, intent(in) :: eps !! Coordinate convergence tolerance; default 0.001.
      integer, optional, intent(in) :: max_iter !! Maximum coordinate sweeps over the path; default 10000.
      real(dp), optional, intent(in) :: group_multiplier(:) !! Optional penalty multiplier per positive group.
      real(dp), allocatable :: lam(:), xw(:,:), yw(:), score(:)
      type(grpreg_preprocess_type) :: prep
      integer :: nl, i, j, g
      real(dp) :: lmin, a, gam, dlt, lmax, factor, z
      character(len=16) :: fam
      fam = 'gaussian'
      if (present(family)) fam = trim(family)
      nl = 100
      if (present(nlambda)) nl = nlambda
      lmin = merge(0.001_dp,0.05_dp,size(x,1) > size(x,2))
      if (present(lambda_min)) lmin = lambda_min
      a = 1.0_dp
      if (present(alpha)) a = alpha
      gam = 0.5_dp
      if (present(gamma)) gam = gamma
      dlt = 1.0e-7_dp
      if (present(delta)) dlt = delta
      if (dlt <= 0.0_dp) error stop 'gBridge: delta must be positive'
      if (gam <= 0.0_dp .or. gam >= 1.0_dp) error stop 'gBridge: gamma must lie in (0,1)'
      if (present(lambda)) then
         allocate(lam(size(lambda)))
         lam = lambda
      else
         call prepare_design(x,group,.true.,xw,prep,group_multiplier)
         allocate(yw(size(y)))
         if (trim(fam) == 'gaussian') then
            yw = y-sum(y)/real(size(y),dp)
            factor = 0.35_dp
         else
            yw = y
            factor = 0.20_dp
         end if
         call null_score_glm(xw,yw,prep%group_reduced,fam,score)
         lmax = 0.0_dp
         do j = 1, size(xw,2)
            g = prep%group_reduced(j)
            if (g <= 0) cycle
            z = abs(dot_product(xw(:,j),score)/real(size(y),dp))/max(prep%group_multiplier(g),tiny(1.0_dp))
            lmax = max(lmax,z)
         end do
         lmax = max(lmax,1.0e-12_dp)*factor**(1.0_dp-gam)/(gam*a)
         if (present(lambda_max)) lmax = lambda_max
         allocate(lam(nl))
         if (nl == 1) then
            lam(1) = lmax
         else if (lmin <= tiny(1.0_dp)) then
            do i = 1, nl-1
               lam(i) = exp(log(lmax)+real(i-1,dp)/real(nl-2,dp)*(log(0.001_dp*lmax)-log(lmax)))
            end do
            lam(nl) = 0.0_dp
            lam = lam(nl:1:-1)
         else
            do i = 1, nl
               lam(i) = exp(log(lmax)+real(i-1,dp)/real(nl-1,dp)*(log(lmin*lmax)-log(lmax)))
            end do
            lam = lam(nl:1:-1)
         end if
      end if
      call grpreg_fit(x,y,group,fit,family=fam,penalty='gBridge',lambda=lam,alpha=a,gamma=gam,eps=eps, &
         max_iter=max_iter,group_multiplier=group_multiplier,bridge_delta=dlt)
   end subroutine gbridge_fit

   subroutine grpsurv_fit(x, time, event, group, fit, penalty, lambda, nlambda, lambda_min, alpha, gamma, tau, eps, &
      max_iter, group_multiplier)
      real(dp), intent(in) :: x(:,:) !! Raw survival-model design matrix, observations by predictors.
      real(dp), intent(in) :: time(:) !! Follow-up or event times for each row of x.
      real(dp), intent(in) :: event(:) !! Event indicators, zero for censoring and one for failure.
      integer, intent(in) :: group(:) !! Predictor group labels; zero identifies unpenalized covariates.
      type(grpreg_fit_type), intent(out) :: fit !! Fitted grouped Cox regularization path.
      character(len=*), optional, intent(in) :: penalty !! Penalty name; default grLasso.
      real(dp), optional, intent(in) :: lambda(:) !! Optional user-specified regularization sequence.
      integer, optional, intent(in) :: nlambda !! Number of generated lambda values; default 100.
      real(dp), optional, intent(in) :: lambda_min !! Minimum generated lambda fraction; default based on n and p.
      real(dp), optional, intent(in) :: alpha !! Sparsity-versus-ridge mixing proportion; default one.
      real(dp), optional, intent(in) :: gamma !! MCP/SCAD concavity parameter.
      real(dp), optional, intent(in) :: tau !! Group exponential lasso tuning parameter; default one third.
      real(dp), optional, intent(in) :: eps !! Coordinate convergence tolerance; default 0.001.
      integer, optional, intent(in) :: max_iter !! Maximum coordinate sweeps over the path; default 10000.
      real(dp), optional, intent(in) :: group_multiplier(:) !! Optional penalty multiplier per positive group.
      integer, allocatable :: ord(:), gnorm(:)
      real(dp), allocatable :: xs(:,:), ts(:), ds(:), xw(:,:), score(:), lam(:), bw(:,:), dev(:), dfr(:), etaw(:,:)
      real(dp), allocatable :: beta(:,:), intercept0(:), intercept_out(:)
      integer, allocatable :: it(:)
      type(grpreg_preprocess_type) :: prep
      character(len=16) :: pen
      real(dp) :: a, gam, tv, tol, lmin
      integer :: nl, mit, ng, i
      logical :: bilevel
      if (size(x,1) /= size(time) .or. size(time) /= size(event)) error stop 'grpsurv: incompatible input lengths'
      pen = 'grLasso'
      if (present(penalty)) pen = trim(penalty)
      if (.not. valid_penalty(pen,.false.)) error stop 'grpsurv: unsupported penalty'
      a = 1.0_dp
      if (present(alpha)) a = alpha
      gam = merge(4.0_dp,3.0_dp,trim(pen) == 'grSCAD')
      if (present(gamma)) gam = gamma
      tv = 1.0_dp/3.0_dp
      if (present(tau)) tv = tau
      tol = 0.001_dp
      if (present(eps)) tol = eps
      mit = 10000
      if (present(max_iter)) mit = max_iter
      nl = 100
      if (present(nlambda)) nl = nlambda
      lmin = merge(0.001_dp,0.05_dp,size(x,1) > size(x,2))
      if (present(lambda_min)) lmin = lambda_min
      call order_real(time,ord)
      allocate(xs(size(x,1),size(x,2)),ts(size(time)),ds(size(event)))
      do i = 1, size(ord)
         xs(i,:) = x(ord(i),:)
         ts(i) = time(ord(i))
         ds(i) = event(ord(i))
      end do
      bilevel = pen(1:2) /= 'gr'
      call prepare_design(xs,group,bilevel,xw,prep,group_multiplier)
      if (present(lambda)) then
         allocate(lam(size(lambda)))
         lam = lambda
      else
         call cox_null_score(ds,score)
         call make_lambda_sequence(xw,score,prep%group_reduced,prep%group_multiplier,a,bilevel,nl,lmin,.true.,lam)
      end if
      call fit_path_cox(xw,ds,prep%group_reduced,prep%group_multiplier,pen,lam,a,gam,tv,tol,mit,bw,dev,dfr,it,etaw)
      allocate(intercept0(size(lam)))
      intercept0 = 0.0_dp
      call restore_coefficients(bw,intercept0,prep,beta,intercept_out)
      call normalize_groups(group,gnorm,ng)
      fit%family = 'cox'
      fit%penalty = pen
      fit%n = size(x,1)
      fit%p = size(x,2)
      fit%ngroups = ng
      fit%nlambda = size(lam)
      fit%alpha = a
      fit%gamma = gam
      fit%tau = tv
      fit%lambda = lam
      allocate(fit%intercept(size(lam)))
      fit%intercept = intercept_out
      fit%beta = beta
      fit%deviance = dev
      fit%df = dfr
      fit%iter = it
      allocate(fit%eta(size(x,1),size(lam)))
      fit%eta = matmul(xs,beta)
      do i = 1, size(lam)
         fit%eta(:,i) = fit%eta(:,i)-sum(fit%eta(:,i))/real(size(x,1),dp)
      end do
      fit%group = gnorm
      fit%group_multiplier = prep%group_multiplier
      fit%x_center = prep%center
      fit%x_scale = prep%scale
      fit%time = ts
      fit%fail = ds
   end subroutine grpsurv_fit

   subroutine coef_grpreg(fit, lambda, beta, intercept)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted grpreg or grpsurv path.
      real(dp), intent(in) :: lambda !! Requested lambda within the fitted path range.
      real(dp), allocatable, intent(out) :: beta(:) !! Interpolated predictor coefficients at lambda.
      real(dp), intent(out) :: intercept !! Interpolated intercept; zero for Cox models.
      integer :: l, r
      real(dp) :: w
      call interpolation_indices(fit%lambda,lambda,l,r,w)
      allocate(beta(fit%p))
      beta = (1.0_dp-w)*fit%beta(:,l)+w*fit%beta(:,r)
      intercept = (1.0_dp-w)*fit%intercept(l)+w*fit%intercept(r)
   end subroutine coef_grpreg

   subroutine predict_grpreg(fit, x, lambda, prediction, type)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted non-survival grpreg path.
      real(dp), intent(in) :: x(:,:) !! Raw predictor matrix at which predictions are required.
      real(dp), intent(in) :: lambda !! Requested lambda within the fitted path range.
      real(dp), allocatable, intent(out) :: prediction(:) !! Link, response, or binary class predictions.
      character(len=*), optional, intent(in) :: type !! Prediction type: link, response, or class; default link.
      real(dp), allocatable :: b(:)
      real(dp) :: b0
      character(len=16) :: typ
      typ = 'link'
      if (present(type)) typ = trim(type)
      call coef_grpreg(fit,lambda,b,b0)
      allocate(prediction(size(x,1)))
      prediction = b0+matmul(x,b)
      if (trim(typ) == 'link' .or. trim(fit%family) == 'gaussian') return
      if (trim(fit%family) == 'binomial') then
         if (trim(typ) == 'response') then
            prediction = logistic(prediction)
         else if (trim(typ) == 'class') then
            where (prediction > 0.0_dp)
               prediction = 1.0_dp
            elsewhere
               prediction = 0.0_dp
            end where
         end if
      else if (trim(fit%family) == 'poisson') then
         if (trim(typ) == 'response') prediction = exp(min(prediction,40.0_dp))
      end if
   end subroutine predict_grpreg

   subroutine predict_grpsurv_link(fit, x, lambda, eta)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted grpsurv Cox path.
      real(dp), intent(in) :: x(:,:) !! Raw predictor matrix for new survival subjects.
      real(dp), intent(in) :: lambda !! Requested lambda within the fitted path range.
      real(dp), allocatable, intent(out) :: eta(:) !! Cox linear predictors for the new subjects.
      real(dp), allocatable :: b(:)
      real(dp) :: b0
      call coef_grpreg(fit,lambda,b,b0)
      allocate(eta(size(x,1)))
      eta = matmul(x,b)
   end subroutine predict_grpsurv_link

   subroutine breslow_baseline(fit, lambda, event_time, cumulative_hazard)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted grpsurv object retaining ordered times and failures.
      real(dp), intent(in) :: lambda !! Requested lambda within the fitted path range.
      real(dp), allocatable, intent(out) :: event_time(:) !! Unique failure times in ascending order.
      real(dp), allocatable, intent(out) :: cumulative_hazard(:) !! Breslow cumulative baseline hazard at event_time.
      integer :: l, r, i, k, ne
      real(dp) :: w, risk, cum
      real(dp), allocatable :: et(:), tmp_t(:), tmp_h(:)
      call interpolation_indices(fit%lambda,lambda,l,r,w)
      allocate(et(fit%n),tmp_t(fit%n),tmp_h(fit%n))
      et = (1.0_dp-w)*fit%eta(:,l)+w*fit%eta(:,r)
      ne = 0
      cum = 0.0_dp
      i = 1
      do while (i <= fit%n)
         if (fit%fail(i) > 0.5_dp) then
            risk = sum(exp(min(et(i:),40.0_dp)))
            k = i
            do while (k <= fit%n .and. abs(fit%time(k)-fit%time(i)) <= 1.0e-12_dp)
               k = k+1
            end do
            cum = cum+sum(fit%fail(i:k-1))/max(risk,tiny(1.0_dp))
            ne = ne+1
            tmp_t(ne) = fit%time(i)
            tmp_h(ne) = cum
            i = k
         else
            i = i+1
         end if
      end do
      allocate(event_time(ne),cumulative_hazard(ne))
      event_time = tmp_t(1:ne)
      cumulative_hazard = tmp_h(1:ne)
   end subroutine breslow_baseline

   subroutine predict_grpsurv_survival(fit, x, lambda, times, survival)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted grpsurv Cox path containing ordered training times and predictors.
      real(dp), intent(in) :: x(:,:) !! Raw predictor matrix for new survival subjects.
      real(dp), intent(in) :: lambda !! Requested lambda within the fitted path range.
      real(dp), intent(in) :: times(:) !! Times at which survival probabilities are requested.
      real(dp), allocatable, intent(out) :: survival(:,:) !! Survival probabilities, subjects by requested times.
      real(dp), allocatable :: et(:), bt(:), bs(:), bh(:)
      integer :: i, j, k
      call predict_grpsurv_link(fit,x,lambda,et)
      call kp_baseline(fit,lambda,bt,bs,bh)
      allocate(survival(size(x,1),size(times)))
      do j = 1, size(times)
         k = step_index(bt,times(j))
         do i = 1, size(x,1)
            survival(i,j) = bs(k)**exp(min(et(i),40.0_dp))
         end do
      end do
   end subroutine predict_grpsurv_survival

   subroutine predict_grpsurv_hazard(fit, x, lambda, times, hazard)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted grpsurv Cox path containing ordered training times and predictors.
      real(dp), intent(in) :: x(:,:) !! Raw predictor matrix for new survival subjects.
      real(dp), intent(in) :: lambda !! Requested lambda within the fitted path range.
      real(dp), intent(in) :: times(:) !! Times at which cumulative hazards are requested.
      real(dp), allocatable, intent(out) :: hazard(:,:) !! Cumulative hazards, subjects by requested times.
      real(dp), allocatable :: et(:), bt(:), bs(:), bh(:)
      integer :: i, j, k
      call predict_grpsurv_link(fit,x,lambda,et)
      call kp_baseline(fit,lambda,bt,bs,bh)
      allocate(hazard(size(x,1),size(times)))
      do j = 1, size(times)
         k = step_index(bt,times(j))
         do i = 1, size(x,1)
            hazard(i,j) = bh(k)*exp(min(et(i),40.0_dp))
         end do
      end do
   end subroutine predict_grpsurv_hazard

   subroutine predict_grpsurv_median(fit, x, lambda, median_time)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted grpsurv Cox path containing ordered training times and predictors.
      real(dp), intent(in) :: x(:,:) !! Raw predictor matrix for new survival subjects.
      real(dp), intent(in) :: lambda !! Requested lambda within the fitted path range.
      real(dp), allocatable, intent(out) :: median_time(:) !! First time survival drops below 0.5; huge if absent.
      real(dp), allocatable :: et(:), bt(:), bs(:), bh(:)
      real(dp) :: s
      integer :: i, k
      call predict_grpsurv_link(fit,x,lambda,et)
      call kp_baseline(fit,lambda,bt,bs,bh)
      allocate(median_time(size(x,1)))
      median_time = huge(1.0_dp)
      do i = 1, size(x,1)
         do k = 2, size(bt)
            s = bs(k)**exp(min(et(i),40.0_dp))
            if (s < 0.5_dp) then
               median_time(i) = bt(k)
               exit
            end if
         end do
      end do
   end subroutine predict_grpsurv_median

   subroutine kp_baseline(fit, lambda, times, survival0, hazard0)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted Cox path used to estimate the upstream Kalbfleisch-Prentice baseline.
      real(dp), intent(in) :: lambda !! Requested lambda at which the baseline is interpolated.
      real(dp), allocatable, intent(out) :: times(:) !! Step-function coordinates beginning at zero then ordered training times.
      real(dp), allocatable, intent(out) :: survival0(:) !! Baseline survival step values at times.
      real(dp), allocatable, intent(out) :: hazard0(:) !! Upstream cumulative-hazard approximation at times.
      real(dp), allocatable :: w(:), risk(:), a(:)
      real(dp) :: frac, base
      integer :: left, right, i
      call interpolation_indices(fit%lambda,lambda,left,right,frac)
      allocate(w(fit%n),risk(fit%n),a(fit%n),times(fit%n+1),survival0(fit%n+1),hazard0(fit%n+1))
      w = (1.0_dp-frac)*exp(min(fit%eta(:,left),40.0_dp))+frac*exp(min(fit%eta(:,right),40.0_dp))
      risk(fit%n) = w(fit%n)
      do i = fit%n-1, 1, -1
         risk(i) = risk(i+1)+w(i)
      end do
      a = 1.0_dp
      do i = 1, fit%n
         if (fit%fail(i) > 0.5_dp .and. risk(i) > w(i)) then
            base = max(1.0_dp-w(i)/risk(i),0.0_dp)
            a(i) = base**(1.0_dp/max(w(i),tiny(1.0_dp)))
         else if (fit%fail(i) > 0.5_dp .and. risk(i) <= w(i)) then
            a(i) = 0.0_dp
         end if
      end do
      times(1) = 0.0_dp
      times(2:) = fit%time
      survival0(1) = 1.0_dp
      hazard0(1) = 0.0_dp
      do i = 1, fit%n
         survival0(i+1) = survival0(i)*a(i)
         hazard0(i+1) = hazard0(i)+(1.0_dp-a(i))
      end do
   end subroutine kp_baseline

   pure integer function step_index(times, query) result(index)
      real(dp), intent(in) :: times(:) !! Ordered step-function coordinates beginning at zero.
      real(dp), intent(in) :: query !! Query time whose constant-step index is requested.
      integer :: i
      index = 1
      do i = 2, size(times)
         if (times(i) <= query) index = i
      end do
   end function step_index

   subroutine loglik_grpreg(fit, loglik, df, active, reml)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted grpreg or grpsurv path.
      real(dp), allocatable, intent(out) :: loglik(:) !! Log likelihood at every fitted lambda.
      real(dp), allocatable, intent(out) :: df(:) !! Effective parameter count associated with loglik.
      logical, optional, intent(in) :: active !! True counts active coefficients instead of stored effective df.
      logical, optional, intent(in) :: reml !! Gaussian-only restricted scale estimate switch.
      logical :: use_active, use_reml
      real(dp) :: rdf, rss
      integer :: l
      use_active = .false.
      if (present(active)) use_active = active
      use_reml = .false.
      if (present(reml)) use_reml = reml
      allocate(loglik(fit%nlambda),df(fit%nlambda))
      if (use_active) then
         do l = 1, fit%nlambda
            df(l) = real(count(abs(fit%beta(:,l)) > 1.0e-12_dp),dp)+merge(0.0_dp,1.0_dp,trim(fit%family)=='cox')
         end do
      else
         df = fit%df
      end if
      if (trim(fit%family) == 'gaussian') then
         do l = 1, fit%nlambda
            rdf = merge(real(fit%n,dp)-df(l),real(fit%n,dp),use_reml)
            rss = max(fit%deviance(l),tiny(1.0_dp))
            loglik(l) = -0.5_dp*real(fit%n,dp)*(log(2.0_dp*acos(-1.0_dp))+log(rss)-log(max(rdf,1.0_dp))) &
               -0.5_dp*rdf
            df(l) = df(l)+1.0_dp
         end do
      else if (trim(fit%family) == 'poisson') then
         do l = 1, fit%nlambda
            loglik(l) = poisson_loglik(fit%y,fit%eta(:,l))
         end do
      else
         loglik = -0.5_dp*fit%deviance
      end if
   end subroutine loglik_grpreg

   subroutine residuals_grpreg(fit, residuals)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted grpreg path retaining responses and linear predictors.
      real(dp), allocatable, intent(out) :: residuals(:,:) !! Gaussian raw or GLM deviance residuals by lambda.
      real(dp) :: mu, d
      integer :: i, l
      allocate(residuals(fit%n,fit%nlambda))
      if (trim(fit%family) == 'cox') then
         call cox_deviance_residuals(fit,residuals)
      else if (trim(fit%family) == 'gaussian') then
         do l = 1, fit%nlambda
            residuals(:,l) = fit%y-fit%eta(:,l)
         end do
      else if (trim(fit%family) == 'binomial') then
         do l = 1, fit%nlambda
            do i = 1, fit%n
               mu = min(max(logistic(fit%eta(i,l)),1.0e-15_dp),1.0_dp-1.0e-15_dp)
               if (fit%y(i) > 0.5_dp) then
                  d = -2.0_dp*log(mu)
               else
                  d = -2.0_dp*log(1.0_dp-mu)
               end if
               residuals(i,l) = sign(sqrt(max(d,0.0_dp)),fit%y(i)-mu)
            end do
         end do
      else
         do l = 1, fit%nlambda
            do i = 1, fit%n
               mu = exp(min(fit%eta(i,l),40.0_dp))
               if (fit%y(i) > 0.0_dp) then
                  d = 2.0_dp*(fit%y(i)*log(fit%y(i)/mu)-(fit%y(i)-mu))
               else
                  d = 2.0_dp*mu
               end if
               residuals(i,l) = sign(sqrt(max(d,0.0_dp)),fit%y(i)-mu)
            end do
         end do
      end if
   end subroutine residuals_grpreg

   subroutine select_grpreg(fit, result, criterion, active)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted regularization path to score by an information criterion.
      type(grpreg_selection_type), intent(out) :: result !! Selected path location, coefficients, df, and criterion vector.
      character(len=*), optional, intent(in) :: criterion !! BIC, AIC, GCV, AICc, or EBIC; default BIC.
      logical, optional, intent(in) :: active !! True uses active coefficient counts rather than stored effective df.
      real(dp), allocatable :: ll(:), dfr(:)
      real(dp) :: pcount, jcount
      character(len=8) :: crit
      integer :: l
      crit = 'BIC'
      if (present(criterion)) crit = trim(criterion)
      call loglik_grpreg(fit,ll,dfr,active=active)
      allocate(result%criterion(fit%nlambda))
      do l = 1, fit%nlambda
         select case (trim(crit))
         case ('AIC')
            result%criterion(l) = -2.0_dp*ll(l)+2.0_dp*dfr(l)
         case ('GCV')
            result%criterion(l) = (-2.0_dp*ll(l)/real(fit%n,dp))/ &
               max((1.0_dp-dfr(l)/real(fit%n,dp))**2,tiny(1.0_dp))
         case ('AICc')
            result%criterion(l) = -2.0_dp*ll(l)+2.0_dp*dfr(l)+ &
               2.0_dp*dfr(l)*(dfr(l)+1.0_dp)/max(real(fit%n,dp)-dfr(l)-1.0_dp,1.0_dp)
         case ('EBIC')
            pcount = real(fit%p,dp)
            jcount = max(0.0_dp,min(pcount,dfr(l)-merge(2.0_dp,1.0_dp,trim(fit%family)=='gaussian')))
            result%criterion(l) = -2.0_dp*ll(l)+log(real(fit%n,dp))*dfr(l)+ &
               2.0_dp*(log_gamma(pcount+1.0_dp)-log_gamma(jcount+1.0_dp)-log_gamma(pcount-jcount+1.0_dp))
         case default
            result%criterion(l) = -2.0_dp*ll(l)+log(real(fit%n,dp))*dfr(l)
         end select
      end do
      result%index = minloc(result%criterion,dim=1)
      result%lambda = fit%lambda(result%index)
      result%df = dfr(result%index)
      allocate(result%beta(fit%p+1))
      result%beta(1) = fit%intercept(result%index)
      result%beta(2:) = fit%beta(:,result%index)
   end subroutine select_grpreg

   pure integer function count_nonzero(fit, index) result(value)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted path whose active predictors are counted.
      integer, intent(in) :: index !! One-based lambda index.
      value = count(abs(fit%beta(:,index)) > 1.0e-12_dp)
   end function count_nonzero

   pure integer function count_nonzero_groups(fit, index) result(value)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted path whose active predictor groups are counted.
      integer, intent(in) :: index !! One-based lambda index.
      integer :: g
      value = 0
      do g = 1, fit%ngroups
         if (any(abs(fit%beta(:,index)) > 1.0e-12_dp .and. fit%group == g)) value = value+1
      end do
   end function count_nonzero_groups

   subroutine group_norms(fit, index, norms)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted path whose raw-scale group norms are requested.
      integer, intent(in) :: index !! One-based lambda index.
      real(dp), allocatable, intent(out) :: norms(:) !! L2 norm of coefficients in each positive group.
      integer :: g
      allocate(norms(fit%ngroups))
      do g = 1, fit%ngroups
         norms(g) = sqrt(sum(fit%beta(:,index)**2,mask=fit%group == g))
      end do
   end subroutine group_norms

   subroutine cox_deviance_residuals(fit, residuals)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted Cox path with ordered times, failures, and linear predictors.
      real(dp), intent(out) :: residuals(:,:) !! Cox deviance residual matrix, observations by lambda.
      real(dp), allocatable :: bt(:), bs(:), bh(:)
      real(dp) :: cumulative, martingale, inside
      integer :: l, i, k
      do l = 1, fit%nlambda
         call kp_baseline(fit,fit%lambda(l),bt,bs,bh)
         do i = 1, fit%n
            k = step_index(bt,fit%time(i))
            cumulative = bh(k)*exp(min(fit%eta(i,l),40.0_dp))
            martingale = fit%fail(i)-cumulative
            if (fit%fail(i) > 0.5_dp) then
               inside = -2.0_dp*(martingale+log(max(fit%fail(i)-martingale,tiny(1.0_dp))))
            else
               inside = -2.0_dp*martingale
            end if
            residuals(i,l) = sign(sqrt(max(inside,0.0_dp)),martingale)
            if (cumulative <= tiny(1.0_dp)) residuals(i,l) = 0.0_dp
         end do
      end do
   end subroutine cox_deviance_residuals

   subroutine null_score_glm(x, y, group, family, score)
      real(dp), intent(in) :: x(:,:) !! Working design matrix for null-model residual computation.
      real(dp), intent(in) :: y(:) !! Centered Gaussian response or raw GLM response.
      integer, intent(in) :: group(:) !! Working group labels with zero marking unpenalized columns.
      character(len=*), intent(in) :: family !! Gaussian, binomial, or Poisson family name.
      real(dp), allocatable, intent(out) :: score(:) !! Null-model score vector used for lambda maximum.
      real(dp), allocatable :: b(:), eta(:), mu(:), r(:)
      real(dp) :: b0, shift, old, maxchg, v, ybar
      integer :: j, it
      allocate(b(size(x,2)),eta(size(y)),mu(size(y)),r(size(y)),score(size(y)))
      b = 0.0_dp
      if (trim(family) == 'gaussian') then
         b0 = 0.0_dp
         eta = 0.0_dp
         do it = 1, 1000
            r = y-eta
            maxchg = 0.0_dp
            do j = 1, size(b)
               if (group(j) /= 0) cycle
               old = b(j)
               b(j) = dot_product(x(:,j),r)/real(size(y),dp)+old
               eta = eta+x(:,j)*(b(j)-old)
               maxchg = max(maxchg,abs(b(j)-old))
            end do
            if (maxchg < 1.0e-10_dp) exit
         end do
         score = y-eta
      else
         ybar = sum(y)/real(size(y),dp)
         if (trim(family) == 'binomial') then
            ybar = min(max(ybar,1.0e-8_dp),1.0_dp-1.0e-8_dp)
            b0 = log(ybar/(1.0_dp-ybar))
         else
            b0 = log(max(ybar,1.0e-12_dp))
         end if
         eta = b0
         do it = 1, 1000
            if (trim(family) == 'binomial') then
               mu = logistic(eta)
               v = 0.25_dp
            else
               mu = exp(min(eta,40.0_dp))
               v = max(maxval(mu),1.0e-8_dp)
            end if
            r = (y-mu)/v
            shift = sum(r)/real(size(y),dp)
            b0 = b0+shift
            eta = eta+shift
            maxchg = abs(shift)
            do j = 1, size(b)
               if (group(j) /= 0) cycle
               old = b(j)
               b(j) = dot_product(x(:,j),r)/real(size(y),dp)+old
               shift = b(j)-old
               r = r-x(:,j)*shift
               eta = eta+x(:,j)*shift
               maxchg = max(maxchg,abs(shift))
            end do
            if (maxchg < 1.0e-10_dp) exit
         end do
         if (trim(family) == 'binomial') then
            score = y-logistic(eta)
         else
            score = y-exp(min(eta,40.0_dp))
         end if
      end if
   end subroutine null_score_glm

   subroutine cox_null_score(event, score)
      real(dp), intent(in) :: event(:) !! Event indicators ordered by increasing follow-up time.
      real(dp), allocatable, intent(out) :: score(:) !! Cox null-model score residuals.
      real(dp) :: cumulative
      integer :: i, n
      n = size(event)
      allocate(score(n))
      cumulative = 0.0_dp
      do i = 1, n
         cumulative = cumulative+event(i)/real(n-i+1,dp)
         score(i) = event(i)-cumulative
      end do
   end subroutine cox_null_score

   subroutine interpolation_indices(path, lambda, left, right, weight)
      real(dp), intent(in) :: path(:) !! Fitted lambda path, normally monotone descending or ascending.
      real(dp), intent(in) :: lambda !! Requested lambda within the path range.
      integer, intent(out) :: left !! First bracketing path index.
      integer, intent(out) :: right !! Second bracketing path index.
      real(dp), intent(out) :: weight !! Linear interpolation weight on the right index.
      integer :: i
      if (lambda < minval(path)-1.0e-12_dp .or. lambda > maxval(path)+1.0e-12_dp) &
         error stop 'grpreg: requested lambda lies outside fitted path'
      left = 1
      right = 1
      weight = 0.0_dp
      do i = 1, size(path)-1
         if ((lambda <= path(i) .and. lambda >= path(i+1)) .or. &
             (lambda >= path(i) .and. lambda <= path(i+1))) then
            left = i
            right = i+1
            if (abs(path(right)-path(left)) > tiny(1.0_dp)) &
               weight = (lambda-path(left))/(path(right)-path(left))
            return
         end if
      end do
      if (abs(lambda-path(size(path))) <= 1.0e-12_dp) then
         left = size(path)
         right = left
      end if
   end subroutine interpolation_indices

   subroutine order_real(x, order)
      real(dp), intent(in) :: x(:) !! Real values to order increasingly.
      integer, allocatable, intent(out) :: order(:) !! Stable one-based permutation producing increasing x.
      integer :: i, j, key
      allocate(order(size(x)))
      order = [(i,i=1,size(x))]
      do i = 2, size(x)
         key = order(i)
         j = i-1
         do while (j >= 1)
            if (x(order(j)) <= x(key)) exit
            order(j+1) = order(j)
            j = j-1
         end do
         order(j+1) = key
      end do
   end subroutine order_real

   pure logical function valid_family(family) result(ok)
      character(len=*), intent(in) :: family !! Candidate grpreg family name.
      ok = trim(family) == 'gaussian' .or. trim(family) == 'binomial' .or. trim(family) == 'poisson'
   end function valid_family

   pure logical function valid_penalty(penalty, allow_bridge) result(ok)
      character(len=*), intent(in) :: penalty !! Candidate grouped or bi-level penalty name.
      logical, intent(in) :: allow_bridge !! Whether standalone group bridge is allowed in this call.
      ok = trim(penalty) == 'grLasso' .or. trim(penalty) == 'grMCP' .or. trim(penalty) == 'grSCAD' .or. &
         trim(penalty) == 'gel' .or. trim(penalty) == 'cMCP'
      if (allow_bridge) ok = ok .or. trim(penalty) == 'gBridge'
   end function valid_penalty

   pure real(dp) function poisson_loglik(y, eta) result(value)
      real(dp), intent(in) :: y(:) !! Nonnegative Poisson outcomes.
      real(dp), intent(in) :: eta(:) !! Log-mean predictors corresponding to y.
      integer :: i
      value = 0.0_dp
      do i = 1, size(y)
         value = value+y(i)*eta(i)-exp(min(eta(i),40.0_dp))-log_gamma(y(i)+1.0_dp)
      end do
   end function poisson_loglik

end module grpreg_api
