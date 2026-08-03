! SPDX-License-Identifier: GPL-3.0-or-later
program test_cla_reference
   use kind_mod, only: dp
   use cla, only: cla_result_t, critical_line
   implicit none

   real(dp), parameter :: tol = 2.0e-10_dp
   real(dp) :: mu(3), covar(3,3), lower(3), upper(3)
   real(dp) :: expected_weights(3,4), expected_lambda(4), expected_sigma(4), expected_mu(4)
   type(cla_result_t) :: result
   integer :: j

   mu = [0.0408_dp, 0.102_dp, -0.023_dp]
   covar = reshape([0.00648_dp,0.00792_dp,0.00473_dp, &
                    0.00792_dp,0.0334_dp,0.0121_dp, &
                    0.00473_dp,0.0121_dp,0.0793_dp],[3,3])
   lower = 0.0_dp
   upper = 1.0_dp
   expected_weights = reshape([ &
      0.0_dp,1.0_dp,0.0_dp, &
      0.991971183148553_dp,0.008028816851447_dp,0.0_dp, &
      0.996984027705433_dp,0.0_dp,0.003015972294567_dp, &
      0.977070230607966_dp,0.0_dp,0.022929769392034_dp],[3,4])
   expected_lambda = [0.416339869281046_dp,0.026683214985438_dp, &
      0.023821645681483_dp,0.0_dp]
   expected_sigma = [0.182756668824971_dp,0.080651550863012_dp, &
      0.080437169955640_dp,0.080248818705100_dp]
   expected_mu = [0.102_dp,0.041291363591309_dp, &
      0.040607580967607_dp,0.039337080712788_dp]

   result = critical_line(mu,covar,lower,upper)
   call assert_true(result%info == 0,'CLA returned failure')
   call assert_true(result%n_turning == 4,'unexpected turning-point count')
   call assert_close(maxval(abs(result%weights-expected_weights)),0.0_dp,tol,'weights')
   call assert_close(maxval(abs(result%lambdas-expected_lambda)),0.0_dp,tol,'lambdas')
   call assert_close(maxval(abs(result%sigma-expected_sigma)),0.0_dp,tol,'sigma')
   call assert_close(maxval(abs(result%mu-expected_mu)),0.0_dp,tol,'mu')
   do j = 1, result%n_turning
      call assert_close(sum(result%weights(:,j)),1.0_dp,1.0e-12_dp,'budget')
      call assert_true(all(result%weights(:,j) >= -1.0e-12_dp),'lower bound')
      call assert_true(all(result%weights(:,j) <= 1.0_dp+1.0e-12_dp),'upper bound')
   end do
   call assert_true(all(result%lambdas(2:) <= result%lambdas(:3)+1.0e-12_dp), &
      'lambdas must decrease')
   print '(a)', 'test_cla_reference: PASS'

contains

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual-expected) > tolerance) then
         write(*,'(a,2es24.14)') trim(label)//' failed: ',actual,expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*,'(a)') trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true

end program test_cla_reference
