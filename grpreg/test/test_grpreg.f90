program test_grpreg
   use grpreg
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   integer, parameter :: n = 40, p = 4
   real(dp) :: x(n,p), yg(n), yb(n), yp(n), time(n), event(n), lam(3)
   integer :: group(p), i, j
   type(grpreg_fit_type) :: fit, fit2
   type(grpreg_cv_type) :: cv
   type(grpreg_selection_type) :: sel
   type(mfdr_result_type) :: mfd
   type(spline_expansion_type) :: spl, nspl
   real(dp), allocatable :: pred(:), res(:,:), ll(:), dfr(:), beta(:), expanded(:,:), gx(:,:), gy(:), gmu(:)
   real(dp), allocatable :: gx2(:,:), gy2(:), gmu2(:), surv(:,:), hazard(:,:), med(:), auc(:), rsq(:), snr(:), sigma(:)
   integer, allocatable :: nv(:), ng(:)
   real(dp) :: b0
   character(len=8), parameter :: penalties(5) = [character(len=8) :: 'grLasso','grMCP','grSCAD','gel','cMCP']

   do i = 1, n
      x(i,1) = (real(i,dp)-20.5_dp)/12.0_dp
      x(i,2) = sin(0.37_dp*real(i,dp))
      x(i,3) = cos(0.23_dp*real(i,dp))
      x(i,4) = sin(0.11_dp*real(i*i,dp))
   end do
   group = [1,1,2,2]
   yg = 1.25_dp+1.7_dp*x(:,1)-0.9_dp*x(:,2)+0.6_dp*x(:,3)-0.3_dp*x(:,4)
   lam = [0.15_dp,0.03_dp,1.0e-8_dp]
   call grpreg_fit(x,yg,group,fit,family='gaussian',penalty='grLasso',lambda=lam,eps=1.0e-9_dp)
   call predict_grpreg(fit,x,lam(3),pred,type='response')
   call assert_true(maxval(abs(pred-yg)) < 5.0e-5_dp,'Gaussian near-unpenalized fit')
   call assert_true(all(fit%deviance >= 0.0_dp),'Gaussian deviance nonnegative')
   call loglik_grpreg(fit,ll,dfr)
   call assert_true(all(ieee_is_finite(ll)),'Gaussian log likelihood finite')
   call residuals_grpreg(fit,res)
   call assert_close(sum(res(:,3)**2),fit%deviance(3),1.0e-8_dp,'Gaussian residual RSS')
   call select_grpreg(fit,sel,criterion='BIC')
   call assert_true(sel%index >= 1 .and. sel%index <= fit%nlambda,'BIC selection index')
   call mfdr_grpreg(fit,mfd)
   call assert_true(all(mfd%mfdr >= 0.0_dp .and. mfd%mfdr <= 1.0_dp),'Gaussian mFDR range')

   do j = 1, size(penalties)
      call grpreg_fit(x,yg,group,fit2,family='gaussian',penalty=trim(penalties(j)),lambda=lam(1:2),eps=1.0e-7_dp)
      call assert_true(all(ieee_is_finite(fit2%deviance)),'Penalty path finite: '//trim(penalties(j)))
   end do
   call gbridge_fit(x,yg,group,fit2,family='gaussian',lambda=[0.01_dp,0.05_dp],eps=1.0e-7_dp)
   call assert_true(trim(fit2%penalty) == 'gBridge','gBridge family wrapper')
   call gbridge_fit(x,yg,group,fit2,family='gaussian',nlambda=5,eps=1.0e-7_dp)
   call assert_true(fit2%lambda(1) <= fit2%lambda(5),'gBridge automatic path is ascending')

   do i = 1, n
      yb(i) = merge(1.0_dp,0.0_dp,mod(i,5) == 0 .or. mod(i,7) == 0 .or. x(i,2) > 0.65_dp)
      yp(i) = real(mod(i,4),dp)
   end do
   call grpreg_fit(x,yb,group,fit2,family='binomial',penalty='grLasso',lambda=[0.05_dp,0.01_dp])
   call predict_grpreg(fit2,x,0.01_dp,pred,type='response')
   call assert_true(all(pred > 0.0_dp .and. pred < 1.0_dp),'Binomial response predictions')
   call mfdr_grpreg(fit2,mfd,x)
   call assert_true(all(mfd%mfdr >= 0.0_dp .and. mfd%mfdr <= 1.0_dp),'Binomial mFDR range')
   call grpreg_fit(x,yp,group,fit2,family='poisson',penalty='grMCP',lambda=[0.05_dp,0.01_dp])
   call predict_grpreg(fit2,x,0.01_dp,pred,type='response')
   call assert_true(all(pred > 0.0_dp),'Poisson response predictions')
   call cv_grpreg(x,yb,group,cv,family='binomial',penalty='grLasso',nfolds=4,nlambda=3)
   call assert_true(size(cv%prediction_error) == 3,'Binomial CV classification metric')
   call cv_grpreg(x,yp,group,cv,family='poisson',penalty='grLasso',nfolds=4,nlambda=3)
   call assert_true(all(cv%cve >= 0.0_dp),'Poisson CV deviance nonnegative')

   call cv_grpreg(x,yg,group,cv,family='gaussian',penalty='grLasso',nfolds=5,nlambda=6)
   call assert_true(size(cv%cve) == 6,'Gaussian CV path length')
   call assert_true(cv%min_index >= 1 .and. cv%min_index <= 6,'Gaussian CV selection')
   call coef_cv_grpreg(cv,beta,b0)
   call assert_true(size(beta) == p,'CV coefficient size')
   call summarize_cv_grpreg(cv,rsq,snr,sigma,nv,ng)
   call assert_true(all(rsq >= 0.0_dp .and. rsq <= 1.0_dp),'CV R-squared range')

   do i = 1, n
      time(i) = real(mod(17*i,41)+1,dp)+0.01_dp*real(i,dp)
      event(i) = merge(1.0_dp,0.0_dp,mod(5*i+2,7) < 4)
   end do
   call grpsurv_fit(x,time,event,group,fit2,penalty='grLasso',lambda=[0.04_dp,0.01_dp])
   call assert_true(trim(fit2%family) == 'cox','Cox family tag')
   call assert_true(all(fit2%deviance >= 0.0_dp),'Cox deviance nonnegative')
   call predict_grpsurv_survival(fit2,x(1:2,:),0.01_dp,[0.0_dp,10.0_dp,30.0_dp],surv)
   call assert_true(all(surv >= 0.0_dp .and. surv <= 1.0_dp),'Survival probability range')
   call assert_true(all(surv(:,2) <= surv(:,1)+1.0e-12_dp),'Survival monotonic first step')
   call predict_grpsurv_hazard(fit2,x(1:2,:),0.01_dp,[0.0_dp,10.0_dp,30.0_dp],hazard)
   call assert_true(all(hazard >= 0.0_dp),'Cumulative hazard nonnegative')
   call predict_grpsurv_median(fit2,x(1:2,:),0.01_dp,med)
   call residuals_grpreg(fit2,res)
   call assert_true(all(ieee_is_finite(res)),'Cox deviance residuals finite')
   call mfdr_grpreg(fit2,mfd,x)
   call assert_true(all(mfd%mfdr >= 0.0_dp .and. mfd%mfdr <= 1.0_dp),'Cox mFDR range')
   call cv_grpsurv(x,time,event,group,cv,penalty='grLasso',nfolds=4,nlambda=4)
   call auc_cv_grpsurv(cv,auc)
   call assert_true(all(auc >= 0.0_dp .and. auc <= 1.0_dp),'Cox CV AUC range')

   call expand_spline(x(:,1:2),spl,df=4,degree=3,spline_type='bs')
   call assert_true(size(spl%x,2) == 8,'B-spline expansion dimensions')
   call predict_spline(spl,x(:,1:2),expanded)
   call assert_true(maxval(abs(expanded-spl%x)) < 1.0e-12_dp,'B-spline prediction reproduces training basis')
   call expand_spline(x(:,1:2),nspl,df=4,spline_type='ns')
   call predict_spline(nspl,x(:,1:2),expanded)
   call assert_true(maxval(abs(expanded-nspl%x)) < 1.0e-9_dp,'Natural-spline prediction reproduces training basis')

   call gen_nonlinear_data(25,8,gx,gy,gmu,seed=1234)
   call gen_nonlinear_data(25,8,gx2,gy2,gmu2,seed=1234)
   call assert_true(maxval(abs(gx-gx2)) <= tiny(1.0_dp),'Generator deterministic seed X')
   call assert_true(maxval(abs(gy-gy2)) <= tiny(1.0_dp),'Generator deterministic seed y')
   call assert_true(maxval(abs(gmu-gmu2)) <= tiny(1.0_dp),'Generator deterministic seed mean')

   print '(a)', 'All grpreg deterministic tests passed.'

contains

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Boolean assertion condition that must hold.
      character(len=*), intent(in) :: label !! Human-readable deterministic-test label.
      if (.not. condition) then
         write(*,'(a)') 'FAILED: '//trim(label)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Computed scalar value being checked.
      real(dp), intent(in) :: expected !! Reference scalar value for the deterministic check.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute difference.
      character(len=*), intent(in) :: label !! Human-readable deterministic-test label.
      if (abs(actual-expected) > tolerance) then
         write(*,'(a,2es20.10)') 'FAILED: '//trim(label)//' actual/reference=',actual,expected
         error stop 1
      end if
   end subroutine assert_close

end program test_grpreg
