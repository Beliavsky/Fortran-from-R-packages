program test_univariate
   use urca_kinds, only : dp
   use urca_types, only : ur_test_result
   use urca_unitroot
   implicit none
   integer, parameter :: n = 240
   real(dp) :: d(n,4), y(n)
   integer :: i, iu
   type(ur_test_result) :: r

   open(newunit=iu,file='test/reference_data.txt',status='old',action='read')
   do i=1,n
      read(iu,*) d(i,:)
   end do
   close(iu)
   y=d(:,1)

   r=adf_test(y,UR_DRIFT,3,LAG_FIXED)
   call assert_zero(r%info,'ADF status')
   call assert_close(r%statistic(1),-4.6794330253858885_dp,2.0e-11_dp,'ADF tau')
   call assert_close(r%statistic(2),10.997867559051583_dp,2.0e-10_dp,'ADF phi1')
   call assert_int(r%lags,3,'ADF lag')

   r=kpss_test(y,KPSS_MU,1)
   call assert_zero(r%info,'KPSS status')
   call assert_close(r%statistic(1),0.0922438642502086_dp,2.0e-12_dp,'KPSS mu')

   r=pp_test(y,PP_ZTAU,PP_CONSTANT,1)
   call assert_zero(r%info,'PP status')
   call assert_close(r%statistic(1),-6.311792252486218_dp,2.0e-11_dp,'PP Z-tau')
   call assert_close(r%auxiliary_statistics(1),-1.8712757523965402_dp,2.0e-11_dp,'PP Z-tau-mu')

   r=ers_test(y,ERS_DFGLS,PP_CONSTANT,4)
   call assert_zero(r%info,'ERS constant DF-GLS status')
   call assert_close(r%statistic(1),-4.220338063518298_dp,2.0e-11_dp,'ERS constant DF-GLS')
   call assert_int(r%lags,4,'ERS constant lag')

   r=ers_test(y,ERS_DFGLS,PP_TREND,4)
   call assert_zero(r%info,'ERS trend DF-GLS status')
   call assert_close(r%statistic(1),-4.277327382260179_dp,2.0e-11_dp,'ERS trend DF-GLS')

   r=ers_test(y,ERS_PTEST,PP_CONSTANT,4)
   call assert_zero(r%info,'ERS P status')
   call assert_close(r%statistic(1),0.6484089160489469_dp,2.0e-11_dp,'ERS P')
   call assert_int(r%lags,1,'ERS P selected lag')

   r=schmidt_phillips_test(y,SP_TAU,3)
   call assert_zero(r%info,'SP tau status')
   call assert_close(r%statistic(1),-6.591188634767892_dp,5.0e-8_dp,'SP tau')
   call assert_int(r%lags,14,'SP HAC lag')

   r=schmidt_phillips_test(y,SP_RHO,3)
   call assert_zero(r%info,'SP rho status')
   call assert_close(r%statistic(1),-75.3611543977774_dp,5.0e-7_dp,'SP rho')

   r=zivot_andrews_test(y,ZA_BOTH,2)
   call assert_zero(r%info,'ZA status')
   call assert_close(r%statistic(1),-5.25698543661746_dp,2.0e-10_dp,'ZA statistic')
   call assert_int(r%break_point,105,'ZA break')

   print '(a)', 'test_univariate: PASS'
contains
   subroutine assert_close(actual,expected,tol,name)
      real(dp), intent(in) :: actual, expected, tol
      character(len=*), intent(in) :: name
      if(abs(actual-expected)>tol)then
         write(*,'(a,2(1x,es24.16))') 'FAIL '//trim(name),actual,expected
         error stop 1
      end if
   end subroutine
   subroutine assert_zero(actual,name)
      integer, intent(in) :: actual
      character(len=*), intent(in) :: name
      if(actual/=0)then
         write(*,'(a,1x,i0)') 'FAIL '//trim(name),actual
         error stop 1
      end if
   end subroutine
   subroutine assert_int(actual,expected,name)
      integer, intent(in) :: actual, expected
      character(len=*), intent(in) :: name
      if(actual/=expected)then
         write(*,'(a,2(1x,i0))') 'FAIL '//trim(name),actual,expected
         error stop 1
      end if
   end subroutine
end program test_univariate
