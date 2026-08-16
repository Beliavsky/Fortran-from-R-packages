program test_models
   use vgam
   implicit none
   integer, parameter :: n=30
   real(dp) :: x(n,2), y(n), t(n), ym(n,2)
   integer :: yc(18), yo(n)
   real(dp), allocatable :: cmat(:,:), pred(:), pmat(:,:)
   logical :: parallel(2)
   type(vglm_result_t) :: gf, pf
   type(multinomial_result_t) :: mf
   type(ordinal_result_t) :: of
   type(beta_regression_result_t) :: bf
   type(negative_binomial_result_t) :: nf
   type(zip_result_t) :: zf
   type(constrained_vglm_result_t) :: cf
   type(vgam_smooth_result_t) :: sf
   integer :: i, failures

   failures=0
   do i=1,n
      t(i)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
      x(i,:)=[1.0_dp,t(i)]
      y(i)=2.0_dp+3.0_dp*t(i)
   end do
   call fit_gaussian(y,x,gf)
   call check(gf%converged,"gaussian convergence",failures)
   call check_close(gf%coefficients(1),2.0_dp,1.0e-9_dp,"gaussian intercept",failures)
   call check_close(gf%coefficients(2),3.0_dp,1.0e-9_dp,"gaussian slope",failures)

   y=exp(0.2_dp+0.5_dp*t)
   call fit_poisson(y,x,pf)
   call check(pf%converged,"poisson convergence",failures)
   call check_close(pf%coefficients(1),0.2_dp,5.0e-7_dp,"poisson intercept",failures)
   call check_close(pf%coefficients(2),0.5_dp,5.0e-7_dp,"poisson slope",failures)

   do i=1,18
      yc(i)=mod(i-1,3)+1
   end do
   call fit_multinomial(yc,x(1:18,:),3,mf)
   call check(mf%converged,"multinomial convergence",failures)
   call check(maxval(abs(sum(mf%fitted_probabilities,dim=2)-1.0_dp))<1.0e-12_dp, &
      "multinomial probabilities",failures)

   do i=1,n
      yo(i)=mod(i-1,3)+1
   end do
   call fit_ordinal(yo,x,3,of)
   call check(of%cutpoints(2)>of%cutpoints(1),"ordinal cutpoint order",failures)
   call check(maxval(abs(sum(of%fitted_probabilities,dim=2)-1.0_dp))<1.0e-10_dp, &
      "ordinal probabilities",failures)

   y=0.15_dp+0.70_dp/(1.0_dp+exp(-(0.4_dp+0.8_dp*t))) + 0.015_dp*sin(7.0_dp*t)
   call fit_beta_regression(y,x,bf)
   call check(bf%precision>0.0_dp,"beta precision",failures)
   call check(all(bf%fitted>0.0_dp.and.bf%fitted<1.0_dp),"beta fitted range",failures)

   do i=1,n
      y(i)=real(mod(3*i+mod(i,4),7),dp)
   end do
   call fit_negative_binomial(y,x,nf)
   call check(nf%size>0.0_dp,"negative binomial size",failures)
   call check(all(nf%fitted>0.0_dp),"negative binomial fitted",failures)

   do i=1,n
      if(mod(i,4)==0)then
         yo(i)=0
      else
         yo(i)=mod(i,3)
      end if
   end do
   call fit_zero_inflated_poisson(yo,x,x(:,1:1),zf,max_iter=200)
   call check(all(zf%zero_probability>0.0_dp.and.zf%zero_probability<1.0_dp), &
      "ZIP zero probabilities",failures)
   call check(all(zf%fitted_mean>=0.0_dp),"ZIP fitted means",failures)

   do i=1,n
      ym(i,1)=1.0_dp+2.0_dp*t(i)
      ym(i,2)=3.0_dp+2.0_dp*t(i)
   end do
   parallel=[.false.,.true.]
   call parallel_constraint(2,2,parallel,cmat)
   call fit_constrained_vglm(ym,x,[family_gaussian,family_gaussian],cmat,cf)
   call check(cf%converged,"constrained VGLM convergence",failures)
   call check_close(cf%coefficients(2,1),2.0_dp,1.0e-8_dp,"parallel slope 1",failures)
   call check_close(cf%coefficients(2,2),2.0_dp,1.0e-8_dp,"parallel slope 2",failures)

   y=sin(2.0_dp*pi*t)+0.05_dp*cos(11.0_dp*t)
   call fit_gam_gaussian(t,y,sf,df=7,lambda=0.1_dp)
   call check(sf%fit%status==0.or.sf%fit%status==100,"spline GAM status",failures)
   pred=sf%predict(t,response=.true.)
   call check(size(pred)==n,"spline prediction size",failures)
   call check(sum((pred-y)**2)/real(n,dp)<0.15_dp,"spline fit quality",failures)

   if(failures/=0)then
      print '(a,i0)',"test_models failures: ",failures
      error stop 1
   end if
   print '(a)',"test_models: PASS"
contains
   subroutine check(ok,name,failures)
      logical,intent(in)::ok
      character(*),intent(in)::name
      integer,intent(inout)::failures
      if(.not.ok)then
         print '(a)',trim(name)//" FAIL"
         failures=failures+1
      end if
   end subroutine check
   subroutine check_close(actual,expected,tol,name,failures)
      real(dp),intent(in)::actual,expected,tol
      character(*),intent(in)::name
      integer,intent(inout)::failures
      call check(abs(actual-expected)<=tol,name,failures)
      if(abs(actual-expected)>tol)print '(a,2(1x,es16.8))',trim(name),actual,expected
   end subroutine check_close
end program test_models
