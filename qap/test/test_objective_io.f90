program test_objective_io
   use qap
   implicit none
   type(qap_problem_t) :: p
   real(dp) :: z

   call read_qaplib('data/qaplib/had12.dat', p)
   if (size(p%A,1) /= 12) error stop 'wrong dimension'
   if (.not. p%has_solution) error stop 'missing solution'
   if (nint(p%opt) /= 1652) error stop 'wrong optimum'
   z = qap_obj(p%A, p%B, p%solution)
   if (abs(z - p%opt) > 1.0e-10_dp) error stop 'objective mismatch'
   if (.not. qap_is_permutation(p%solution)) error stop 'invalid solution permutation'
   print *, 'test_objective_io: PASS'
end program test_objective_io
