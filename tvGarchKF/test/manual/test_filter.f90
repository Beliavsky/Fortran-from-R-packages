program test_filter
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use fgarch_kinds, only : dp
   use tvgarchkf
   use test_support
   implicit none
   type(tvgarch_spec) :: spec, bad_spec
   type(tvgarch_filter_result) :: result
   real(dp) :: y(3), expected_p(3), expected_f(3), expected_k(3), expected_x(3), expected_v(3)
   real(dp) :: y_missing(3)

   spec = make_tvgarch_spec(make_tv_function([0.1_dp]), &
                            make_tv_function([0.2_dp]), &
                            make_tv_function([0.6_dp]))
   y = [0.0_dp,1.0_dp,-0.5_dp]
   result = tvgarch_kalman_filter(y,spec,corrected_constraints=.true.)
   call assert_true(result%status == 0,'filter status')
   expected_p = [0.11111111111111117_dp,0.036000000000000025_dp,0.012509652509652514_dp]
   expected_f = [1.1111111111111112_dp,1.036_dp,1.0125096525096524_dp]
   expected_k = [0.26_dp,0.22084942084942089_dp,0.20741305674191585_dp]
   expected_x = [0.0_dp,-0.13_dp,0.03513513513513518_dp]
   expected_v = [0.5_dp,0.37_dp,0.5351351351351351_dp]
   call assert_all_close(result%state_variance,expected_p,1.0e-12_dp,'state variance')
   call assert_all_close(result%mse,expected_f,1.0e-12_dp,'mse')
   call assert_all_close(result%gain,expected_k,1.0e-12_dp,'gain')
   call assert_all_close(result%state,expected_x,1.0e-12_dp,'state')
   call assert_all_close(result%conditional_variance,expected_v,1.0e-12_dp,'conditional variance')
   call assert_close(result%criterion,0.831790649644333_dp,1.0e-12_dp,'criterion')

   y_missing = [0.0_dp,ieee_value(0.0_dp,ieee_quiet_nan),0.5_dp]
   result = tvgarch_kalman_filter(y_missing,spec,predict=2,corrected_constraints=.true.)
   call assert_true(result%status == 0,'missing filter status')
   call assert_true(size(result%sigma) == 5,'forecast length')
   call assert_true(all(result%conditional_variance > 0.0_dp),'positive forecast variance')

   bad_spec = make_tvgarch_spec(make_tv_function([-0.1_dp]), &
                                make_tv_function([0.2_dp]), &
                                make_tv_function([0.6_dp]))
   result = tvgarch_kalman_filter(y,bad_spec,corrected_constraints=.true.)
   call assert_true(result%status /= 0,'corrected constraints reject negative omega')
   write(*,'(a)') 'test_filter: PASS'
end program test_filter
