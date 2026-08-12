program solve_had20
   use qap
   implicit none
   type(qap_problem_t) :: p
   type(qap_control_t) :: ctl
   type(qap_result_t) :: res

   call read_qaplib('data/qaplib/had20.dat', p)
   ctl%rep = 10
   ctl%seed = 1000_i64
   ctl%verbose = .true.
   call qap_solve(p%A, p%B, res, ctl)

   write(*,'(/,a,f12.0)') 'best objective: ', res%objective
   write(*,'(a,f12.0)') 'known optimum:  ', p%opt
   write(*,'(a,20(i0,1x))') 'permutation:    ', res%permutation
end program solve_had20
