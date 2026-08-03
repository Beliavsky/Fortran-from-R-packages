! SPDX-License-Identifier: Artistic-2.0
program test_distribution
   use ldhmm
   implicit none
   type(ecld_type) :: distribution
   real(dp) :: dx, integral, x
   integer :: i, status

   distribution = ecld_create(lambda=1.0_dp, sigma=1.0_dp, mu=0.0_dp, status=status)
   call assert_true(status == LDHMM_SUCCESS, 'constructor status')
   call assert_close(ecld_pdf(distribution,0.0_dp), 1.0_dp/sqrt(acos(-1.0_dp)), &
      1.0e-13_dp, 'normal-form density')
   call assert_close(ecld_cdf(distribution,0.0_dp), 0.5_dp, 1.0e-14_dp, 'cdf at mean')
   call assert_close(ecld_variance(distribution), 0.5_dp, 1.0e-13_dp, 'variance')
   call assert_close(ecld_sd(distribution), sqrt(0.5_dp), 1.0e-13_dp, 'sd')
   call assert_close(ecld_skewness(distribution), 0.0_dp, 1.0e-15_dp, 'skewness')
   call assert_close(ecld_kurtosis(distribution), 3.0_dp, 1.0e-12_dp, 'kurtosis')
   call assert_close(ecld_cdf(distribution,1.25_dp)+ecld_cdf(distribution,-1.25_dp), &
      1.0_dp, 1.0e-12_dp, 'cdf symmetry')

   dx = 0.001_dp
   integral = 0.0_dp
   do i = 0, 12000
      x = -6.0_dp + real(i,dp)*dx
      if (i == 0 .or. i == 12000) then
         integral = integral + 0.5_dp*ecld_pdf(distribution,x)*dx
      else
         integral = integral + ecld_pdf(distribution,x)*dx
      end if
   end do
   call assert_close(integral, 1.0_dp, 2.0e-10_dp, 'density integral')

   distribution = ecld_create(lambda=3.0_dp, sigma=0.02_dp, mu=0.001_dp)
   call assert_true(ecld_pdf(distribution,distribution%mu) > 0.0_dp, 'heavy-tail density')
   call assert_close(ecld_ccdf(distribution,0.01_dp), &
      1.0_dp-ecld_cdf(distribution,0.01_dp), 1.0e-14_dp, 'ccdf')

   print '(a)', 'test_distribution: PASS'

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*,'(a)') 'FAIL: '//message
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: message
      call assert_true(abs(actual-expected) <= tolerance, message)
   end subroutine assert_close

end program test_distribution
