program test_sa_had20
   use qap
   implicit none
   type(qap_problem_t) :: p
   type(qap_control_t) :: ctl
   type(qap_result_t) :: res
   real(dp) :: z

   call read_qaplib('data/qaplib/had20.dat', p)
   ctl%rep = 10
   ctl%seed = 1000_i64
   call qap_solve(p%A, p%B, res, ctl)

   if (.not. qap_is_permutation(res%permutation)) error stop 'invalid SA permutation'
   z = qap_obj(p%A, p%B, res%permutation)
   if (abs(z-res%objective) > 1.0e-8_dp) error stop 'stored objective mismatch'
   if (abs(res%objective-p%opt) > 1.0e-8_dp) error stop 'had20 optimum not reached'
   if (res%attempted_swaps <= 0_i64) error stop 'no SA work recorded'
   if (res%accepted_swaps <= 0_i64) error stop 'no accepted swaps recorded'
   print *, 'test_sa_had20: PASS'
end program test_sa_had20
