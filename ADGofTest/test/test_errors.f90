program test_errors
   use adgoftest, only : dp, ad_statistic, ad_test_result, ad_test_uniform, &
      ad_empty_sample, ad_probability_out_of_range, ad_success
   implicit none

   real(dp), allocatable :: empty(:)
   real(dp), parameter :: endpoint_sample(3) = [0.0_dp, 0.5_dp, 1.0_dp]
   real(dp), parameter :: invalid_sample(3) = [-0.1_dp, 0.5_dp, 0.9_dp]
   type(ad_test_result) :: result
   real(dp) :: value
   integer :: status

   allocate(empty(0))
   call ad_test_uniform(empty, result)
   call assert_true(result%status == ad_empty_sample, 'empty sample status')

   value = ad_statistic(invalid_sample, status=status)
   call assert_true(status == ad_probability_out_of_range, 'invalid range status')
   call assert_true(abs(value) <= tiny(1.0_dp), 'invalid range value')

   call ad_test_uniform(endpoint_sample, result)
   call assert_true(result%status == ad_success, 'endpoint status')
   call assert_true(result%statistic > huge(1.0_dp), 'endpoint infinite statistic')
   call assert_true(abs(result%p_value) <= tiny(1.0_dp), 'endpoint p-value')

   call ad_test_uniform(endpoint_sample, result, clip_probabilities=.true., clamp_p_value=.true.)
   call assert_true(result%status == ad_success, 'clipped endpoint status')
   call assert_true(result%statistic > 0.0_dp, 'clipped endpoint statistic')
   call assert_true(result%p_value >= 0.0_dp .and. result%p_value <= 1.0_dp, 'clipped endpoint p-value')

contains

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a)') trim(label)
         error stop 1
      end if
   end subroutine assert_true

end program test_errors
