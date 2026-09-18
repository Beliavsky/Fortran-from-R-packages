module grpreg_cv
   use grpreg_kinds, only : dp
   use grpreg_types, only : grpreg_cv_type, grpreg_fit_type
   use grpreg_api, only : grpreg_fit, grpsurv_fit, predict_grpreg, coef_grpreg, count_nonzero, count_nonzero_groups
   use grpreg_math, only : logistic
   implicit none
   private
   public :: cv_grpreg, cv_grpsurv, coef_cv_grpreg, predict_cv_grpreg, auc_cv_grpsurv, summarize_cv_grpreg

contains

   subroutine cv_grpreg(x, y, group, cv, family, penalty, nfolds, fold, nlambda, lambda_min, alpha, gamma, tau)
      real(dp), intent(in) :: x(:,:) !! Raw predictor matrix, observations by predictors, used for all fold fits.
      real(dp), intent(in) :: y(:) !! Response vector for Gaussian, binomial, or Poisson cross-validation.
      integer, intent(in) :: group(:) !! Predictor group labels; zero marks unpenalized predictors.
      type(grpreg_cv_type), intent(out) :: cv !! Cross-validation result with full-data fit and fold predictions.
      character(len=*), optional, intent(in) :: family !! Outcome family; gaussian is used when omitted.
      character(len=*), optional, intent(in) :: penalty !! Group penalty name; grLasso is used when omitted.
      integer, optional, intent(in) :: nfolds !! Number of deterministic folds when fold is omitted; default ten.
      integer, optional, intent(in) :: fold(:) !! Optional one-based fold assignment for each observation.
      integer, optional, intent(in) :: nlambda !! Number of automatically generated lambda values; default one hundred.
      real(dp), optional, intent(in) :: lambda_min !! Smallest lambda as a fraction of lambda maximum.
      real(dp), optional, intent(in) :: alpha !! Group-penalty mixing fraction in (0,1].
      real(dp), optional, intent(in) :: gamma !! MCP/SCAD concavity parameter.
      real(dp), optional, intent(in) :: tau !! Group exponential-lasso decay parameter.
      type(grpreg_fit_type) :: ffold
      character(len=16) :: fam, pen
      integer, allocatable :: f(:), train(:), test(:)
      real(dp), allocatable :: xtr(:,:), ytr(:), xte(:,:), pred(:)
      real(dp) :: a, gam, tv, lmin, mu, err, mean_e, var_e
      integer :: nf, nl, i, j, k, nt, nv
      fam = 'gaussian'
      if (present(family)) fam = trim(family)
      pen = 'grLasso'
      if (present(penalty)) pen = trim(penalty)
      nf = 10
      if (present(nfolds)) nf = nfolds
      nl = 100
      if (present(nlambda)) nl = nlambda
      a = 1.0_dp
      if (present(alpha)) a = alpha
      gam = merge(4.0_dp,3.0_dp,trim(pen) == 'grSCAD')
      if (present(gamma)) gam = gamma
      tv = 1.0_dp/3.0_dp
      if (present(tau)) tv = tau
      lmin = merge(0.0001_dp,0.05_dp,size(x,1) > size(x,2))
      if (present(lambda_min)) lmin = lambda_min
      call grpreg_fit(x,y,group,cv%fit,family=fam,penalty=pen,nlambda=nl,lambda_min=lmin,alpha=a,gamma=gam,tau=tv)
      allocate(f(size(y)))
      if (present(fold)) then
         if (size(fold) /= size(y)) error stop 'cv.grpreg: fold length does not match response length'
         f = fold
         nf = maxval(f)
      else
         call make_folds(y,fam,nf,f)
      end if
      cv%nfolds = nf
      cv%fold = f
      cv%lambda = cv%fit%lambda
      allocate(cv%cv_prediction(size(y),cv%fit%nlambda))
      cv%cv_prediction = 0.0_dp
      do i = 1, nf
         nt = count(f /= i)
         nv = count(f == i)
         allocate(train(nt),test(nv),xtr(nt,size(x,2)),ytr(nt),xte(nv,size(x,2)))
         train = pack([(j,j=1,size(y))],f /= i)
         test = pack([(j,j=1,size(y))],f == i)
         do j = 1, nt
            xtr(j,:) = x(train(j),:)
            ytr(j) = y(train(j))
         end do
         do j = 1, nv
            xte(j,:) = x(test(j),:)
         end do
         call grpreg_fit(xtr,ytr,group,ffold,family=fam,penalty=pen,lambda=cv%fit%lambda,alpha=a,gamma=gam,tau=tv)
         do k = 1, cv%fit%nlambda
            call predict_grpreg(ffold,xte,cv%fit%lambda(k),pred,type='response')
            cv%cv_prediction(test,k) = pred
         end do
         deallocate(train,test,xtr,ytr,xte)
      end do
      allocate(cv%cve(cv%fit%nlambda),cv%cvse(cv%fit%nlambda))
      if (trim(fam) == 'binomial') allocate(cv%prediction_error(cv%fit%nlambda))
      do k = 1, cv%fit%nlambda
         mean_e = 0.0_dp
         do i = 1, size(y)
            err = point_deviance(y(i),cv%cv_prediction(i,k),fam)
            mean_e = mean_e+err
         end do
         mean_e = mean_e/real(size(y),dp)
         cv%cve(k) = mean_e
         var_e = 0.0_dp
         if (size(y) > 1) then
            do i = 1, size(y)
               err = point_deviance(y(i),cv%cv_prediction(i,k),fam)
               var_e = var_e+(err-mean_e)**2
            end do
            var_e = var_e/real(size(y)-1,dp)
         end if
         cv%cvse(k) = sqrt(max(var_e,0.0_dp)/real(size(y),dp))
         if (trim(fam) == 'binomial') then
            cv%prediction_error(k) = real(count((cv%cv_prediction(:,k) < 0.5_dp) .eqv. (y > 0.5_dp)),dp)/real(size(y),dp)
         end if
      end do
      cv%min_index = minloc(cv%cve,dim=1)
      cv%lambda_min = cv%lambda(cv%min_index)
      if (trim(fam) == 'gaussian') then
         mu = sum(y)/real(size(y),dp)
         cv%null_deviance = sum((y-mu)**2)/real(size(y),dp)
      else if (trim(fam) == 'binomial') then
         mu = min(max(sum(y)/real(size(y),dp),1.0e-10_dp),1.0_dp-1.0e-10_dp)
         cv%null_deviance = sum([(point_deviance(y(i),mu,fam),i=1,size(y))])/real(size(y),dp)
      else
         mu = max(sum(y)/real(size(y),dp),1.0e-10_dp)
         cv%null_deviance = sum([(point_deviance(y(i),mu,fam),i=1,size(y))])/real(size(y),dp)
      end if
   end subroutine cv_grpreg

   subroutine cv_grpsurv(x, time, event, group, cv, penalty, nfolds, fold, nlambda, lambda_min, alpha, gamma, tau)
      real(dp), intent(in) :: x(:,:) !! Raw predictor matrix for Cox cross-validation.
      real(dp), intent(in) :: time(:) !! Nonnegative follow-up times for survival observations.
      real(dp), intent(in) :: event(:) !! Event indicators, with one for failures and zero for censoring.
      integer, intent(in) :: group(:) !! Predictor group labels; zero marks unpenalized predictors.
      type(grpreg_cv_type), intent(out) :: cv !! Cox cross-validation result with out-of-fold linear predictors.
      character(len=*), optional, intent(in) :: penalty !! Group penalty name; grLasso is used when omitted.
      integer, optional, intent(in) :: nfolds !! Number of deterministic folds when fold is omitted; default ten.
      integer, optional, intent(in) :: fold(:) !! Optional fold assignments in the original observation order.
      integer, optional, intent(in) :: nlambda !! Number of automatically generated lambda values; default one hundred.
      real(dp), optional, intent(in) :: lambda_min !! Smallest lambda as a fraction of lambda maximum.
      real(dp), optional, intent(in) :: alpha !! Group-penalty mixing fraction in (0,1].
      real(dp), optional, intent(in) :: gamma !! MCP/SCAD concavity parameter.
      real(dp), optional, intent(in) :: tau !! Group exponential-lasso decay parameter.
      type(grpreg_fit_type) :: ffold
      character(len=16) :: pen
      integer, allocatable :: f(:), train(:), test(:)
      real(dp), allocatable :: xtr(:,:), ttr(:), dtr(:), xte(:,:), eta(:), loss(:)
      real(dp) :: a, gam, tv, lmin, mean_l, var_l
      integer :: nf, nl, i, j, k, nt, nv, ne
      pen = 'grLasso'
      if (present(penalty)) pen = trim(penalty)
      nf = 10
      if (present(nfolds)) nf = nfolds
      nl = 100
      if (present(nlambda)) nl = nlambda
      a = 1.0_dp
      if (present(alpha)) a = alpha
      gam = merge(4.0_dp,3.0_dp,trim(pen) == 'grSCAD')
      if (present(gamma)) gam = gamma
      tv = 1.0_dp/3.0_dp
      if (present(tau)) tv = tau
      lmin = merge(0.001_dp,0.05_dp,size(x,1) > size(x,2))
      if (present(lambda_min)) lmin = lambda_min
      call grpsurv_fit(x,time,event,group,cv%fit,penalty=pen,nlambda=nl,lambda_min=lmin,alpha=a,gamma=gam,tau=tv)
      allocate(f(size(time)))
      if (present(fold)) then
         if (size(fold) /= size(time)) error stop 'cv.grpsurv: fold length does not match survival data'
         f = fold
         nf = maxval(f)
      else
         call make_folds(event,'binomial',nf,f)
      end if
      cv%nfolds = nf
      cv%fold = f
      cv%lambda = cv%fit%lambda
      allocate(cv%cv_prediction(size(time),cv%fit%nlambda))
      cv%cv_prediction = 0.0_dp
      do i = 1, nf
         nt = count(f /= i)
         nv = count(f == i)
         allocate(train(nt),test(nv),xtr(nt,size(x,2)),ttr(nt),dtr(nt),xte(nv,size(x,2)))
         train = pack([(j,j=1,size(time))],f /= i)
         test = pack([(j,j=1,size(time))],f == i)
         do j = 1, nt
            xtr(j,:) = x(train(j),:)
            ttr(j) = time(train(j))
            dtr(j) = event(train(j))
         end do
         do j = 1, nv
            xte(j,:) = x(test(j),:)
         end do
         call grpsurv_fit(xtr,ttr,dtr,group,ffold,penalty=pen,lambda=cv%fit%lambda,alpha=a,gamma=gam,tau=tv)
         do k = 1, cv%fit%nlambda
            call predict_cox_at(ffold,xte,cv%fit%lambda(k),eta)
            cv%cv_prediction(test,k) = eta
         end do
         deallocate(train,test,xtr,ttr,dtr,xte)
      end do
      allocate(cv%cve(cv%fit%nlambda),cv%cvse(cv%fit%nlambda))
      ne = max(1,count(event > 0.5_dp))
      do k = 1, cv%fit%nlambda
         call cox_event_losses(time,event,cv%cv_prediction(:,k),loss)
         cv%cve(k) = sum(loss)/real(ne,dp)
         mean_l = sum(loss)/real(max(1,size(loss)),dp)
         var_l = 0.0_dp
         if (size(loss) > 1) var_l = sum((loss-mean_l)**2)/real(size(loss)-1,dp)
         cv%cvse(k) = sqrt(max(var_l,0.0_dp)*real(size(loss),dp))/real(ne,dp)
      end do
      cv%min_index = minloc(cv%cve,dim=1)
      cv%lambda_min = cv%lambda(cv%min_index)
      cv%null_deviance = cv%cve(1)
   end subroutine cv_grpsurv

   subroutine coef_cv_grpreg(cv, beta, intercept)
      type(grpreg_cv_type), intent(in) :: cv !! Cross-validation object whose minimum-CV coefficients are requested.
      real(dp), allocatable, intent(out) :: beta(:) !! Predictor coefficients at the minimum-CV lambda.
      real(dp), intent(out) :: intercept !! Intercept at the minimum-CV lambda; zero for Cox fits.
      call coef_grpreg(cv%fit,cv%lambda_min,beta,intercept)
   end subroutine coef_cv_grpreg

   subroutine predict_cv_grpreg(cv, x, prediction, type)
      type(grpreg_cv_type), intent(in) :: cv !! Cross-validation object whose minimum-CV fit is used for prediction.
      real(dp), intent(in) :: x(:,:) !! Raw predictor matrix at which minimum-CV predictions are requested.
      real(dp), allocatable, intent(out) :: prediction(:) !! Predictions at the selected minimum-CV lambda.
      character(len=*), optional, intent(in) :: type !! Link, response, or class prediction type; default link.
      character(len=16) :: typ
      typ = 'link'
      if (present(type)) typ = trim(type)
      if (trim(cv%fit%family) == 'cox') then
         call predict_cox_at(cv%fit,x,cv%lambda_min,prediction)
      else
         call predict_grpreg(cv%fit,x,cv%lambda_min,prediction,type=typ)
      end if
   end subroutine predict_cv_grpreg

   subroutine auc_cv_grpsurv(cv, auc)
      type(grpreg_cv_type), intent(in) :: cv !! Cox cross-validation result containing out-of-fold linear predictors.
      real(dp), allocatable, intent(out) :: auc(:) !! Harrell-style concordance statistic at each fitted lambda.
      integer :: l, i, j
      real(dp) :: good, total, si, sj
      if (trim(cv%fit%family) /= 'cox') error stop 'AUC.cv.grpsurv: fit is not a Cox model'
      allocate(auc(cv%fit%nlambda))
      do l = 1, cv%fit%nlambda
         good = 0.0_dp
         total = 0.0_dp
         do i = 1, cv%fit%n-1
            do j = i+1, cv%fit%n
               if (cv%fit%time(i) < cv%fit%time(j) .and. cv%fit%fail(i) > 0.5_dp) then
                  si = cv%cv_prediction(i,l)
                  sj = cv%cv_prediction(j,l)
               else if (cv%fit%time(j) < cv%fit%time(i) .and. cv%fit%fail(j) > 0.5_dp) then
                  si = cv%cv_prediction(j,l)
                  sj = cv%cv_prediction(i,l)
               else
                  cycle
               end if
               total = total+1.0_dp
               if (si > sj) good = good+1.0_dp
               if (abs(si-sj) <= 1.0e-12_dp) good = good+0.5_dp
            end do
         end do
         if (total > 0.0_dp) then
            auc(l) = good/total
         else
            auc(l) = 0.5_dp
         end if
      end do
   end subroutine auc_cv_grpsurv

   subroutine summarize_cv_grpreg(cv, r_squared, snr, sigma, nvars, ngroups)
      type(grpreg_cv_type), intent(in) :: cv !! Cross-validation object summarized along its lambda sequence.
      real(dp), allocatable, intent(out) :: r_squared(:) !! Cross-validated R-squared or deviance-explained measure.
      real(dp), allocatable, intent(out) :: snr(:) !! Cross-validated signal-to-noise ratio estimate.
      real(dp), allocatable, intent(out) :: sigma(:) !! Gaussian residual-scale estimates; zeros for other families.
      integer, allocatable, intent(out) :: nvars(:) !! Number of active coefficients at each lambda.
      integer, allocatable, intent(out) :: ngroups(:) !! Number of active predictor groups at each lambda.
      real(dp), allocatable :: signal(:)
      integer :: l
      allocate(r_squared(size(cv%cve)),snr(size(cv%cve)),sigma(size(cv%cve)),signal(size(cv%cve)))
      allocate(nvars(size(cv%cve)),ngroups(size(cv%cve)))
      signal = max(cv%null_deviance-cv%cve,0.0_dp)
      if (trim(cv%fit%family) == 'gaussian') then
         r_squared = min(max(1.0_dp-cv%cve/max(cv%null_deviance,tiny(1.0_dp)),0.0_dp),1.0_dp)
         sigma = sqrt(max(cv%cve,0.0_dp))
      else
         r_squared = min(max(1.0_dp-exp(cv%cve-cv%null_deviance),0.0_dp),1.0_dp)
         sigma = 0.0_dp
      end if
      snr = signal/max(cv%cve,tiny(1.0_dp))
      do l = 1, size(cv%cve)
         nvars(l) = count_nonzero(cv%fit,l)
         ngroups(l) = count_nonzero_groups(cv%fit,l)
      end do
   end subroutine summarize_cv_grpreg

   subroutine make_folds(y, family, nfolds, fold)
      real(dp), intent(in) :: y(:) !! Response or censoring indicator used for deterministic fold balancing.
      character(len=*), intent(in) :: family !! Family name; binomial requests class-balanced folds.
      integer, intent(in) :: nfolds !! Positive number of folds to create.
      integer, intent(out) :: fold(:) !! One-based deterministic fold assignment for each observation.
      integer :: i, n0, n1
      n0 = 0
      n1 = 0
      if (trim(family) == 'binomial') then
         do i = 1, size(y)
            if (y(i) > 0.5_dp) then
               n1 = n1+1
               fold(i) = mod(n1-1,nfolds)+1
            else
               n0 = n0+1
               fold(i) = mod(n0+n1-1,nfolds)+1
            end if
         end do
      else
         do i = 1, size(y)
            fold(i) = mod(i-1,nfolds)+1
         end do
      end if
   end subroutine make_folds

   pure real(dp) function point_deviance(y, prediction, family) result(value)
      real(dp), intent(in) :: y !! Observed response value.
      real(dp), intent(in) :: prediction !! Response-scale prediction at the same observation.
      character(len=*), intent(in) :: family !! Gaussian, binomial, or Poisson family name.
      real(dp) :: p, yly
      select case (trim(family))
      case ('gaussian')
         value = (y-prediction)**2
      case ('binomial')
         p = min(max(prediction,1.0e-5_dp),1.0_dp-1.0e-5_dp)
         if (y > 0.5_dp) then
            value = -2.0_dp*log(p)
         else
            value = -2.0_dp*log(1.0_dp-p)
         end if
      case default
         p = max(prediction,1.0e-12_dp)
         yly = 0.0_dp
         if (y > 0.0_dp) yly = y*log(y)
         value = 2.0_dp*(yly-y+p-y*log(p))
      end select
   end function point_deviance

   subroutine predict_cox_at(fit, x, lambda, eta)
      type(grpreg_fit_type), intent(in) :: fit !! Fitted Cox path used for new linear predictors.
      real(dp), intent(in) :: x(:,:) !! Raw predictor matrix for new Cox subjects.
      real(dp), intent(in) :: lambda !! Lambda value at which coefficients are interpolated.
      real(dp), allocatable, intent(out) :: eta(:) !! Interpolated Cox linear predictors.
      real(dp), allocatable :: b(:)
      real(dp) :: b0
      call coef_grpreg(fit,lambda,b,b0)
      allocate(eta(size(x,1)))
      eta = matmul(x,b)
   end subroutine predict_cox_at

   subroutine cox_event_losses(time, event, eta, loss)
      real(dp), intent(in) :: time(:) !! Follow-up times corresponding to eta.
      real(dp), intent(in) :: event(:) !! Event indicators corresponding to eta.
      real(dp), intent(in) :: eta(:) !! Cross-validated Cox linear predictors.
      real(dp), allocatable, intent(out) :: loss(:) !! Per-event negative twice partial-log-likelihood contributions.
      integer, allocatable :: ord(:)
      real(dp), allocatable :: es(:), ts(:), ds(:), risk(:)
      integer :: i, j, ne
      call order_values(time,ord)
      allocate(es(size(eta)),ts(size(time)),ds(size(event)),risk(size(eta)))
      do i = 1, size(time)
         ts(i) = time(ord(i))
         ds(i) = event(ord(i))
         es(i) = eta(ord(i))
      end do
      risk(size(eta)) = exp(min(es(size(eta)),40.0_dp))
      do i = size(eta)-1, 1, -1
         risk(i) = risk(i+1)+exp(min(es(i),40.0_dp))
      end do
      ne = count(ds > 0.5_dp)
      allocate(loss(ne))
      j = 0
      do i = 1, size(eta)
         if (ds(i) > 0.5_dp) then
            j = j+1
            loss(j) = -2.0_dp*(es(i)-log(max(risk(i),tiny(1.0_dp))))
         end if
      end do
   end subroutine cox_event_losses

   subroutine order_values(x, order)
      real(dp), intent(in) :: x(:) !! Real values to sort increasingly.
      integer, allocatable, intent(out) :: order(:) !! One-based stable permutation producing increasing values.
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
   end subroutine order_values

end module grpreg_cv
