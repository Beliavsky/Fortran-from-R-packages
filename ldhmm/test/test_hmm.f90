! SPDX-License-Identifier: Artistic-2.0
program test_hmm
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use ldhmm
   implicit none
   type(ldhmm_model) :: model, decoded
   real(dp) :: param(2,3), gamma_matrix(2,2), delta(2), x(8)
   real(dp), allocatable :: log_alpha(:, :), log_beta(:, :), forecast_state(:, :)
   real(dp), allocatable :: forecast_density(:, :), conditional(:, :), grid(:)
   real(dp), allocatable :: residuals(:), history(:, :)
   real(dp) :: mllk, log_likelihood, dx
   integer, allocatable :: states(:)
   integer :: i, status

   param(1,:) = [0.002936740_dp,0.01977561_dp,1.141693_dp]
   param(2,:) = [-0.001707031_dp,0.03718047_dp,1.324177_dp]
   gamma_matrix(1,:) = [0.98083875_dp,0.01916125_dp]
   gamma_matrix(2,:) = [0.04931245_dp,0.95068755_dp]
   delta = [0.7201662_dp,0.2798338_dp]
   x = [0.01_dp,-0.02_dp,0.005_dp,0.03_dp,-0.01_dp,0.0_dp,0.04_dp,-0.03_dp]
   model = ldhmm_create(2,param,gamma_matrix,delta,status=status)

   mllk = ldhmm_mllk(model,x,status)
   call ldhmm_log_forward(model,x,log_alpha,status)
   call ldhmm_log_backward(model,x,log_beta,status)
   log_likelihood = local_log_sum_exp(log_alpha(:,size(x)))
   call assert_close(mllk,-log_likelihood,1.0e-12_dp,'forward likelihood')
   call assert_true(maxval(abs(log_beta(:,size(x)))) <= 1.0e-15_dp,'backward terminal')
   do i = 1, size(x)
      call assert_close(local_log_sum_exp(log_alpha(:,i)+log_beta(:,i)), &
         log_likelihood,1.0e-11_dp,'forward backward identity')
   end do

   decoded = ldhmm_decode(model,x,status=status)
   do i = 1, size(x)
      call assert_close(sum(decoded%states_prob(:,i)),1.0_dp,1.0e-13_dp, &
         'posterior normalization')
   end do
   states = ldhmm_viterbi(model,x,status)
   call assert_true(size(states) == size(x),'viterbi size')
   call assert_true(all(states >= 1 .and. states <= model%m),'viterbi range')

   forecast_state = ldhmm_forecast_state(model,x,10,status)
   do i = 1, 10
      call assert_close(sum(forecast_state(:,i)),1.0_dp,1.0e-13_dp, &
         'forecast state normalization')
   end do

   allocate(grid(601))
   dx = 0.001_dp
   grid = [(-0.3_dp+dx*real(i-1,dp),i=1,size(grid))]
   conditional = ldhmm_conditional_prob(model,x,grid,status)
   do i = 1, size(x)
      call assert_close(sum(conditional(:,i))*dx,1.0_dp,3.0e-5_dp, &
         'conditional density integral')
   end do
   forecast_density = ldhmm_forecast_prob(model,x,grid,5,status)
   do i = 1, 5
      call assert_close(sum(forecast_density(i,:))*dx,1.0_dp,3.0e-5_dp, &
         'forecast density integral')
   end do

   residuals = ldhmm_pseudo_residuals(model,x,grid_length=300,status=status)
   call assert_true(all(ieee_is_finite(residuals)),'pseudo residuals finite')
   history = ldhmm_decode_stats_history(decoded,moving_average_order=3)
   call assert_true(all(shape(history) == [size(x),3]),'history shape')
   call assert_true(all(ieee_is_finite(history)),'history finite')
   print '(a)', 'test_hmm: PASS'

contains

   real(dp) function local_log_sum_exp(values) result(answer)
      real(dp), intent(in) :: values(:)
      real(dp) :: maximum
      maximum = maxval(values)
      answer = maximum + log(sum(exp(values-maximum)))
   end function local_log_sum_exp

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

end program test_hmm
