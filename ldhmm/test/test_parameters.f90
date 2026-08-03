! SPDX-License-Identifier: Artistic-2.0
program test_parameters
   use ldhmm
   implicit none
   type(ldhmm_model) :: model, recovered, nonstationary, recovered_nonstationary
   real(dp) :: param(3,3), gamma_matrix(3,3), delta(3)
   real(dp), allocatable :: working(:), initialized(:, :), stationary(:)
   integer :: status

   param(1,:) = [0.003_dp,0.020_dp,1.0_dp]
   param(2,:) = [-0.003_dp,0.030_dp,1.1_dp]
   param(3,:) = [-0.006_dp,0.030_dp,1.3_dp]
   gamma_matrix(1,:) = [0.980_dp,0.019_dp,0.001_dp]
   gamma_matrix(2,:) = [0.030_dp,0.960_dp,0.010_dp]
   gamma_matrix(3,:) = [0.001_dp,0.109_dp,0.890_dp]
   delta = [0.596_dp,0.367_dp,0.037_dp]

   model = ldhmm_create(3,param,gamma_matrix,stationary=.true.,status=status)
   call assert_true(status == LDHMM_SUCCESS, 'stationary constructor')
   call assert_close(sum(model%delta),1.0_dp,1.0e-12_dp,'stationary sum')
   call assert_vector_close(matmul(model%delta,model%gamma),model%delta,1.0e-10_dp, &
      'stationary equation')
   working = ldhmm_natural_to_working(model,mu_scale=0.01_dp)
   recovered = ldhmm_working_to_natural(model,working,mu_scale=0.01_dp,status=status)
   call assert_true(status == LDHMM_SUCCESS, 'working transform status')
   call assert_matrix_close(recovered%param,model%param,1.0e-13_dp,'parameter round trip')
   call assert_matrix_close(recovered%gamma,model%gamma,1.0e-13_dp,'gamma round trip')

   nonstationary = ldhmm_create(3,param,gamma_matrix,delta,stationary=.false.,status=status)
   working = ldhmm_natural_to_working(nonstationary)
   recovered_nonstationary = ldhmm_working_to_natural(nonstationary,working,status=status)
   call assert_vector_close(recovered_nonstationary%delta,delta,1.0e-13_dp, &
      'delta round trip')

   initialized = ldhmm_gamma_init(3)
   call assert_vector_close(initialized(1,:),[0.95_dp,0.04_dp,0.01_dp],1.0e-14_dp, &
      'gamma init row 1')
   call assert_vector_close(initialized(2,:),[0.04_dp,0.92_dp,0.04_dp],1.0e-14_dp, &
      'gamma init row 2')
   call assert_vector_close(initialized(3,:),[0.01_dp,0.04_dp,0.95_dp],1.0e-14_dp, &
      'gamma init row 3')

   call ldhmm_stationary_distribution(gamma_matrix,stationary,status)
   call assert_true(status == LDHMM_SUCCESS, 'stationary solver')
   call assert_close(sum(stationary),1.0_dp,1.0e-12_dp,'stationary solver sum')
   print '(a)', 'test_parameters: PASS'

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

   subroutine assert_vector_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      character(len=*), intent(in) :: message
      call assert_true(size(actual) == size(expected), message//' size')
      call assert_true(maxval(abs(actual-expected)) <= tolerance, message)
   end subroutine assert_vector_close

   subroutine assert_matrix_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual(:, :), expected(:, :), tolerance
      character(len=*), intent(in) :: message
      call assert_true(all(shape(actual) == shape(expected)), message//' shape')
      call assert_true(maxval(abs(actual-expected)) <= tolerance, message)
   end subroutine assert_matrix_close

end program test_parameters
