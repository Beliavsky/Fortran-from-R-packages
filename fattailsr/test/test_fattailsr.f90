program test_fattailsr
   use fattailsr
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none
   integer :: failures

   failures = 0
   call test_math(failures)
   call test_parameters(failures)
   call test_distribution(failures)
   call test_tail_risk(failures)
   call test_moments(failures)
   call test_estimation(failures)
   call test_returns(failures)
   call test_random_generation(failures)

   if (failures /= 0) then
      print '(a,i0)', 'FAILED tests: ', failures
      error stop 1
   end if
   print '(a)', 'All FatTailsR tests passed.'

contains

   subroutine assert_close(name, actual, expected, tol, failures)
      character(*), intent(in) :: name
      real(dp), intent(in) :: actual, expected, tol
      integer, intent(inout) :: failures
      real(dp) :: scale
      scale = max(1.0_dp, abs(expected))
      if (.not. ieee_is_finite(actual) .or. abs(actual-expected) > tol*scale) then
         failures = failures + 1
         print '(a,2(1x,es24.15),1x,a,es12.4)', trim(name), actual, expected, 'tol=',tol
      end if
   end subroutine assert_close

   subroutine assert_true(name, condition, failures)
      character(*), intent(in) :: name
      logical, intent(in) :: condition
      integer, intent(inout) :: failures
      if (.not. condition) then
         failures = failures + 1
         print '(a)', trim(name)//' failed'
      end if
   end subroutine assert_true

   subroutine test_math(failures)
      integer, intent(inout) :: failures
      real(dp) :: p
      p = 0.123456789_dp
      call assert_close('invlogit(logit(p))', invlogit(logit(p)), p, 2.0e-14_dp, failures)
      call assert_close('regularized beta', regularized_beta(0.5_dp,2.0_dp,3.0_dp), &
                        0.6875_dp, 1.0e-13_dp, failures)
      call assert_close('normal quantile 0.975', normal_quantile(0.975_dp), &
                        1.959963984540054_dp, 2.0e-9_dp, failures)
      call assert_close('normal cdf inverse', normal_cdf(normal_quantile(0.234_dp)), &
                        0.234_dp, 2.0e-9_dp, failures)
      call assert_close('kashp derivative', dkashp_dx(2.0_dp,3.0_dp), &
                        3.0_dp/sqrt(13.0_dp), 1.0e-14_dp, failures)
   end subroutine test_math

   subroutine test_parameters(failures)
      integer, intent(inout) :: failures
      type(kiener_parameters) :: par
      par = make_k2(0.2_dp,1.1_dp,6.666666666666666_dp,5.454545454545454_dp)
      call assert_close('aw2k',par%k,6.0_dp,1.0e-14_dp,failures)
      call assert_close('aw2e',par%e,0.1_dp,1.0e-14_dp,failures)
      call assert_close('aw2d',par%d,0.1_dp/6.0_dp,1.0e-14_dp,failures)
      call assert_close('kd2a inverse',kd2a(par%k,par%d),par%a,1.0e-14_dp,failures)
      call assert_close('kd2w inverse',kd2w(par%k,par%d),par%w,1.0e-14_dp,failures)
      call assert_close('ke2a inverse',ke2a(par%k,par%e),par%a,1.0e-14_dp,failures)
      call assert_close('kw2e inverse',kw2e(par%k,par%w),par%e,1.0e-14_dp,failures)
      call assert_true('parameters valid',parameters_valid(par),failures)
   end subroutine test_parameters

   subroutine test_distribution(failures)
      integer, intent(inout) :: failures
      type(kiener_parameters) :: par
      real(dp), parameter :: probs(7) = [0.01_dp,0.05_dp,0.25_dp,0.5_dp,0.75_dp,0.95_dp,0.99_dp]
      real(dp), parameter :: expected_q(7) = [ &
         -2.641147808822987_dp,-1.5692416096733304_dp,-0.45783940452728505_dp, &
          0.2_dp,0.8823762862704037_dp,2.1516964166288113_dp,3.5114255956953486_dp]
      real(dp) :: q, p, h
      integer :: i
      par = make_k4(0.2_dp,1.1_dp,6.0_dp,0.1_dp)
      do i=1,size(probs)
         q = qkiener(probs(i),par)
         call assert_close('K4 quantile',q,expected_q(i),3.0e-13_dp,failures)
         p = pkiener(q,par)
         call assert_close('CDF quantile inversion',p,probs(i),2.0e-12_dp,failures)
         call assert_close('density derivative reciprocity',dpkiener(probs(i),par)*&
                           dqkiener(probs(i),par),1.0_dp,2.0e-13_dp,failures)
         call assert_close('x/logit inversion',lkiener(q,par),logit(probs(i)),2.0e-11_dp,failures)
      end do
      call assert_close('K1 closed wrapper',qkiener1(0.9_dp,0.0_dp,1.0_dp,4.0_dp), &
                        qkiener(0.9_dp,make_k1(0.0_dp,1.0_dp,4.0_dp)),1.0e-14_dp,failures)
      h = 1.0e-6_dp
      q = qkiener(0.7_dp,par)
      call assert_close('density finite difference',dkiener(q,par), &
         (pkiener(q+h,par)-pkiener(q-h,par))/(2.0_dp*h),2.0e-7_dp,failures)
      call assert_close('standard logistic sd scale median',qlogisst(0.5_dp,2.0_dp,3.0_dp), &
                        2.0_dp,1.0e-14_dp,failures)
      call assert_close('standard logistic inversion',plogisst(qlogisst(0.13_dp,2.0_dp,3.0_dp), &
                        2.0_dp,3.0_dp),0.13_dp,2.0e-14_dp,failures)
   end subroutine test_distribution

   subroutine test_tail_risk(failures)
      integer, intent(inout) :: failures
      type(kiener_parameters) :: par
      par = make_k4(0.2_dp,1.1_dp,6.0_dp,0.1_dp)
      call assert_close('left tail mean 0.01',ltmkiener(0.01_dp,par), &
                        -3.4062319044517353_dp,2.0e-12_dp,failures)
      call assert_close('right tail mean 0.99',rtmkiener(0.99_dp,par), &
                        4.584749503668668_dp,2.0e-12_dp,failures)
      call assert_close('ES left absolute',eskiener(0.01_dp,par), &
                        3.4062319044517353_dp,2.0e-12_dp,failures)
      call assert_close('ES right',eskiener(0.99_dp,par), &
                        4.584749503668668_dp,2.0e-12_dp,failures)
      call assert_close('signed ES left',eskiener(0.01_dp,par,.true.), &
                        -3.4062319044517353_dp,2.0e-12_dp,failures)
      call assert_true('tail mean below quantile',ltmkiener(0.05_dp,par)<qkiener(0.05_dp,par),failures)
      call assert_true('right tail mean above quantile',rtmkiener(0.95_dp,par)>qkiener(0.95_dp,par),failures)
      call assert_close('c correction median limit',ckiener(0.5_dp,par),1.0_dp,1.0e-13_dp,failures)
      call assert_true('h correction finite',ieee_is_finite(hkiener(0.01_dp,par)),failures)
   end subroutine test_tail_risk

   subroutine test_moments(failures)
      integer, intent(inout) :: failures
      type(kiener_parameters) :: par
      type(moment_summary) :: s
      par = make_k4(0.2_dp,1.1_dp,6.0_dp,0.1_dp)
      call assert_close('raw mean',kiener_raw_moment(1,par), &
                        0.23550060108041992_dp,3.0e-12_dp,failures)
      call assert_close('central variance',kiener_central_moment(2,par), &
                        1.396413049252459_dp,3.0e-12_dp,failures)
      call assert_close('central third',kiener_central_moment(3,par), &
                        0.7033238161164664_dp,4.0e-12_dp,failures)
      call assert_close('central fourth',kiener_central_moment(4,par), &
                        14.4810346978933_dp,5.0e-12_dp,failures)
      s = kiener_moment_summary(par)
      call assert_close('moment sd',s%standard_deviation,1.181699221144052_dp,3.0e-12_dp,failures)
      call assert_close('moment skewness',s%skewness,0.4262206391340166_dp,4.0e-12_dp,failures)
      call assert_close('moment excess',s%excess_kurtosis,4.426288159843972_dp,6.0e-12_dp,failures)
   end subroutine test_moments

   subroutine test_estimation(failures)
      integer, intent(inout) :: failures
      integer, parameter :: n=199
      type(kiener_parameters) :: truth, est5, est7, est11, fitted
      real(dp) :: x(n), p5(5),p7(7),p11(11),x5(5),x7(7),x11(11)
      integer :: i
      truth=make_k4(0.2_dp,1.1_dp,5.0_dp,0.2_dp)
      do i=1,n
         x(i)=qkiener(real(i,dp)/real(n+1,dp),truth)
      end do
      p5=five_probs(n,4)
      do i=1,5; x5(i)=qkiener(p5(i),truth); end do
      est5=estimkiener5(x5,p5,20.0_dp,0.90_dp)
      call assert_close('estim5 m',est5%m,truth%m,1.0e-13_dp,failures)
      call assert_close('estim5 g',est5%g,truth%g,2.0e-5_dp,failures)
      call assert_close('estim5 k',est5%k,truth%k,2.0e-5_dp,failures)
      call assert_close('estim5 e',est5%e,truth%e,2.0e-5_dp,failures)
      p7=seven_probs(n)
      do i=1,7; x7(i)=qkiener(p7(i),truth); end do
      est7=estimkiener7(x7,p7,20.0_dp)
      call assert_close('estim7 m',est7%m,truth%m,1.0e-13_dp,failures)
      call assert_close('estim7 g',est7%g,truth%g,2.0e-2_dp,failures)
      call assert_close('estim7 k',est7%k,truth%k,3.0e-1_dp,failures)
      call assert_close('estim7 e',est7%e,truth%e,3.0e-2_dp,failures)
      p11=eleven_probs(n)
      do i=1,11; x11(i)=qkiener(p11(i),truth); end do
      est11=estimkiener11(x11,p11,7,20.0_dp)
      call assert_close('estim11 m',est11%m,truth%m,1.0e-13_dp,failures)
      call assert_close('estim11 g',est11%g,truth%g,2.0e-2_dp,failures)
      call assert_close('estim11 k',est11%k,truth%k,3.0e-1_dp,failures)
      call assert_close('estim11 e',est11%e,truth%e,3.0e-2_dp,failures)
      fitted=fit_kiener_k4(x,maxk=20.0_dp,mink=1.53_dp,maxe=0.5_dp)
      call assert_close('fit m',fitted%m,truth%m,1.0e-12_dp,failures)
      call assert_close('fit g',fitted%g,truth%g,2.0e-5_dp,failures)
      call assert_close('fit k',fitted%k,truth%k,2.0e-4_dp,failures)
      call assert_close('fit e',fitted%e,truth%e,2.0e-5_dp,failures)
   end subroutine test_estimation

   subroutine test_returns(failures)
      integer, intent(inout) :: failures
      real(dp) :: prices(4),returns(4),missing(4),filled(4),nan
      nan=ieee_value(nan,ieee_quiet_nan)
      prices=[100.0_dp,101.0_dp,100.0_dp,110.0_dp]
      call price_returns(prices,returns)
      call assert_close('return first zero',returns(1),0.0_dp,1.0e-14_dp,failures)
      call assert_close('log return',returns(2),100.0_dp*log(1.01_dp),1.0e-13_dp,failures)
      call price_returns(prices,returns,log_returns=.false.,multiplier=1.0_dp)
      call assert_close('simple return',returns(2),0.01_dp,1.0e-14_dp,failures)
      missing=[nan,2.0_dp,nan,4.0_dp]
      call replace_nonfinite(missing,filled)
      call assert_close('fill leading',filled(1),2.0_dp,1.0e-14_dp,failures)
      call assert_close('fill internal',filled(3),2.0_dp,1.0e-14_dp,failures)
      call assert_close('elevate negative',elevate(-3.0_dp,4.0_dp),1.0_dp,1.0e-14_dp,failures)
   end subroutine test_returns

   subroutine test_random_generation(failures)
      integer, intent(inout) :: failures
      real(dp) :: x(1000),y(1000)
      type(kiener_parameters) :: par
      par=make_k4(0.0_dp,1.0_dp,6.0_dp,0.1_dp)
      call rkiener(x,par)
      call rlogisst(y)
      call assert_true('random Kiener finite',all(ieee_is_finite(x)),failures)
      call assert_true('random logistic finite',all(ieee_is_finite(y)),failures)
      call assert_true('random Kiener nonconstant',maxval(x)>minval(x),failures)
   end subroutine test_random_generation

end program test_fattailsr
