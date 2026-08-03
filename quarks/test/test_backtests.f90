program test_backtests
   use quarks
   implicit none
   real(dp), parameter :: tol = 2.0e-12_dp
   real(dp) :: loss(8), var(8), es(8)
   type(coverage_result) :: coverage
   type(traffic_result) :: traffic
   type(loss_result) :: score

   loss = [1.0_dp, 2.0_dp, 3.0_dp, 1.0_dp, 4.0_dp, 2.0_dp, 5.0_dp, 1.0_dp]
   var = [2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 2.0_dp]
   es = [2.5_dp, 2.5_dp, 2.5_dp, 2.5_dp, 3.5_dp, 3.5_dp, 4.5_dp, 2.5_dp]

   coverage = cvgtest(loss, var, 0.75_dp)
   if (coverage%violations /= 3) error stop 'wrong violation count'
   if (coverage%n00 /= 1 .or. coverage%n01 /= 3 .or. &
       coverage%n10 /= 3 .or. coverage%n11 /= 0) error stop 'wrong transitions'
   call assert_close(coverage%lr_uc, 0.609575080709437_dp, tol, 'LRuc')
   call assert_close(coverage%lr_ind, 6.08633065357725_dp, tol, 'LRind')
   call assert_close(coverage%p_cc, 0.03515625_dp, tol, 'p.cc')

   traffic = trftest(loss, var, 0.75_dp)
   call assert_close(traffic%cumulative_probability, 0.8861846923828125_dp, &
      tol, 'traffic probability')

   score = lossfun(loss, es)
   call assert_close(score%lossfun1, 7500.0_dp, tol, 'loss 1')
   call assert_close(score%lossfun2, 7513.5_dp, tol, 'loss 2')
   call assert_close(score%lossfun3, 7506.5_dp, tol, 'loss 3')
   call assert_close(score%lossfun4, 7506.5_dp, tol, 'loss 4')
   print *, 'test_backtests: PASS'

contains

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual - expected) > tolerance) then
         print *, trim(label), actual, expected
         error stop 1
      end if
   end subroutine assert_close

end program test_backtests
