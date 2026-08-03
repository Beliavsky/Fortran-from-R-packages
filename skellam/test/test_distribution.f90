program test_distribution
   use skellam, only : dp, i8, dskellam, pskellam, skellam_log_pmf
   implicit none

   integer(i8), parameter :: k(6) = [-5_i8, -1_i8, 0_i8, 1_i8, 5_i8, 15_i8]
   real(dp), parameter :: pmf_ref(6) = [ &
      0.004592447845998783_dp, 0.12201556296996496_dp, &
      0.16772188586190170_dp, 0.18302334445494730_dp, &
      0.034873900830553264_dp, 1.0714183521874816e-7_dp]
   real(dp), parameter :: cdf_ref(6) = [ &
      0.006401909316797408_dp, 0.24698869937222823_dp, &
      0.41471058523413000_dp, 0.59773392968907730_dp, &
      0.97555379187455690_dp, 0.99999997627263950_dp]
   real(dp) :: pmf(6), cdf(6), total
   integer :: i
   integer(i8) :: j

   do i = 1, size(k)
      pmf(i) = dskellam(k(i), 3.0_dp, 2.0_dp)
      cdf(i) = pskellam(real(k(i), dp), 3.0_dp, 2.0_dp)
      call assert_close(pmf(i), pmf_ref(i), 2.0e-12_dp, 'PMF reference')
      call assert_close(cdf(i), cdf_ref(i), 2.0e-12_dp, 'CDF reference')
      call assert_close(log(pmf(i)), skellam_log_pmf(k(i), 3.0_dp, 2.0_dp), &
         2.0e-12_dp, 'log PMF')
   end do

   call assert_close(dskellam(3_i8, 5.0_dp, 0.0_dp), &
      exp(-5.0_dp)*5.0_dp**3/6.0_dp, 1.0e-14_dp, 'Poisson positive limit')
   call assert_close(dskellam(-3_i8, 0.0_dp, 5.0_dp), &
      exp(-5.0_dp)*5.0_dp**3/6.0_dp, 1.0e-14_dp, 'Poisson negative limit')
   call assert_close(dskellam(1_i8, 0.0_dp, 5.0_dp), 0.0_dp, 0.0_dp, 'outside support')

   total = 0.0_dp
   do j = -40_i8, 40_i8
      total = total + dskellam(j, 3.0_dp, 2.0_dp)
   end do
   call assert_close(total, 1.0_dp, 2.0e-14_dp, 'PMF normalization')
   print '(a)', 'test_distribution: PASS'

contains

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(*), intent(in) :: label
      if (abs(actual - expected) > tolerance*max(1.0_dp, abs(expected))) then
         write(*,'(a,2(1x,es24.16))') trim(label)//' failed:', actual, expected
         error stop 1
      end if
   end subroutine assert_close

end program test_distribution
